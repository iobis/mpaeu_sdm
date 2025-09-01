############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
#################### Diversity maps - richness and others ######################

library(dplyr)
library(glue)
library(terra)
source("functions/post_div_functions.R")

# Settings
results_folder <- "/data/scps/v5/results"
out_folder <- "/data/scps/v5/diversity"
preproc_folder <- file.path(out_folder, "preproc")
fs::dir_create(out_folder)
fs::dir_create(preproc_folder)
model_acro <- "mpaeu"
target_thresholds <- c("p10", "mss")
target_types <- c("std", "const")

# Step 1 - pre-process layers
done_species <- get_done_species(results_folder)

prepare_layers(
  thresholds = target_thresholds,
  type = target_types,
  results_folder = results_folder,
  out_folder = preproc_folder,
  parallel = TRUE,
  n_cores = 110,
  max_mem = TRUE,
  global_mask = "data/shapefiles/mpa_europe_starea_v3.gpkg",
  base_file = "data/env/current/thetao_baseline_depthsurf_mean.tif",
  species = done_species,
  verbose = TRUE
)

# Step 2 - assemble
species_list <- get_sp_list()
groups <- unique(species_list$group)

proc_list <- list_processed(preproc_folder, write_csv = FALSE)

types_grid <- expand.grid(
  sel_threshold = target_thresholds,
  type = target_types
)

scen_grid <- data.frame(
  scenario = c("current", rep(paste0("ssp", c(126, 245, 370, 460, 585)), 2)),
  decade = c(NA, rep(c(2050, 2100), each = 5))
)

save_raster <- \(r, file, as_cog = TRUE) {
  r <- as.int(r)
  writeRaster(r, file, datatype = "INT2U")
  obissdm::cogeo_optim(file)
  return(invisible(NULL))
}

for (tg in seq_len(nrow(types_grid))) {
  cli::cli_alert_info(cli::bg_cyan("Combination {tg} out of {nrow(types_grid)}"))
  thresh <- types_grid$sel_threshold[tg]
  mtype <- types_grid$type[tg]

  for (g in seq_len(length(groups))) {
    cli::cli_alert_info(cli::bg_green("Group {g} out of {length(groups)}"))
    if (groups[g] == "all") {
      sel_species <- species_list$taxonID
    } else {
      sel_species <- species_list$taxonID[species_list$group == groups[g]]
    }
    av_species <- proc_list |>
      filter(taxonid %in% sel_species) |>
      filter(threshold == thresh) |>
      filter(type == mtype)

    pb <- progress::progress_bar$new(total = nrow(scen_grid))
    for (sc in seq_len(nrow(scen_grid))) {
      pb$tick()
      scen <- scen_grid$scenario[sc]
      decade <- scen_grid$decade[sc]

      if (scenario == "current") {
        scen_list <- av_species |>
          filter(scenario == scen)
      } else {
        scen_list <- av_species |>
          filter(scenario == scen) |>
          filter(period == decade)
      }

      if (nrow(scen_list) != length(unique(scen_list$taxonid))) stop("Potential problem with list. Check.")

      layers <- rast(file.path(preproc_folder, scen_list$file))

      # Produce traditional layer
      cont_richness <- sum(layers)
      cont_richness <- cont_richness / 100

      # Produce binary layer
      bin_layers <- classify(layers, matrix(c(0, Inf, 1), nrow = 1))
      bin_richness <- sum(bin_layers)

      # SESAM probabilities
      sesam_richness <- sesam_prr(layers, cont_richness)
      sesam_richness <- sum(sesam_richness)

      # Save rasters
      base_f <- paste0(
        "metric=richness_",
        gsub("\\.tif", "", gsub("taxonid=.*_model", "model", scen_list$file[1]))
      )
      tempf <- paste0(base_f, "_group=", gsub("\\/", "-", tolower(groups[g])), "_what=continuous.tif")
      save_raster(cont_richness, file.path(out_folder, tempf))
      tempf <- paste0(base_f, "_group=", gsub("\\/", "-", tolower(groups[g])), "_what=binary.tif")
      save_raster(bin_richness, file.path(out_folder, tempf))
      tempf <- paste0(base_f, "_group=", gsub("\\/", "-", tolower(groups[g])), "_what=sesam.tif")
      save_raster(sesam_richness, file.path(out_folder, tempf))
    }

    # Load true richness
    coarse_raster <- aggregate(cont_richness, 10)

    records <- dplyr::bind_rows(
      lapply(
        file.path(
          results_folder,
          glue::glue("taxonid={av_species$taxonid}/model={model_acro}/taxonid={av_species$taxonid}_model={model_acro}_what=fitocc.parquet")
        ),
        \(x) {
          df <- arrow::read_parquet(x)
          df$cell <- cellFromXY(coarse_raster, as.data.frame(df))
          df$taxonID <- gsub("_.*", "", basename(x))
          df <- df[!is.na(df$cell), ]
          df
        }
      )
    )
    records <- records |>
      group_by(cell) |>
      distinct(taxonID) |>
      count() |>
      rename(richness = n)

    true_richness <- coarse_raster
    true_richness[] <- NA
    true_richness[records$cell] <- records$richness
    true_richness <- mask(true_richness, coarse_raster)

    base_f <- paste0(
      "metric=richness_",
      gsub("\\.tif", "", gsub("taxonid=.*_model", "model", scen_list$file[1]))
    )
    tempf <- paste0(base_f, "_group=", gsub("\\/", "-", tolower(groups[g])), "_what=raw.tif")
    save_raster(true_richness, file.path(out_folder, tempf))
  }
}




# OLD CODE
# Load packages/settings -----
library(progress)
library(arrow)
library(dplyr)
library(terra)
library(furrr)
source("functions/post_div_functions.R")
# Beta diversity
# library(DBI)
# library(duckdb)
# library(adespatial)
# library(betapart)


# Settings
batch_size <- 200
acro <- "mpaeu"
results_folder <- "/data/scps/v3/results"
outf <- "/data/scps/v3/diversity"
fs::dir_create(outf)
starea <- "data/shapefiles/mpa_europe_starea_v2.shp"
# Mask target defines which type of mask to use
# 1 = "native_ecoregions"; 2 = "fit_ecoregions"; 3 = "fit_region"; 
# 4 = "convex_hull"; 5 = "minbounding_circle"; 6 = "buffer100m"
mask_target <- 3
# Number of cores for parallel processing
n_cores <- 15
plan(multisession, workers = n_cores)


# List species
cat("\n\nPreparing data...\n")
sp <- list.files(results_folder)

sp_status <- purrr::map(sp, function(id) {
  js <- jsonlite::read_json(file.path(results_folder, glue::glue("{id}/model={acro}/{id}_model={acro}_what=log.json")))
  good_models <- unlist(js$model_good)
  data.frame(
    taxonID = id,
    maxent = ifelse("maxent" %in% good_models, TRUE, FALSE),
    rf = ifelse("rf" %in% good_models, TRUE, FALSE),
    xgboost = ifelse("xgboost" %in% good_models, TRUE, FALSE),
    esm = ifelse(is.numeric(good_models), TRUE, FALSE),
    ensemble = ifelse(length(good_models) > 1, TRUE, FALSE)
  )
}, .progress = T)
sp_status <- do.call("rbind", sp_status)
sp_status$taxonID <- as.numeric(gsub("taxonid=", "", sp_status$taxonID))

priority <- c("ensemble", "maxent", "rf", "xgboost", "esm")

sp_status$which_valid <- apply(sp_status[, priority], 1, function(row) {
  valid_model <- priority[which(row)[1]]
  ifelse(is.na(valid_model), "none", valid_model)
})

groups <- get_sp_list()
groups <- groups[, c("scientificName", "taxonID", "group")]

sp_status <- left_join(sp_status, groups)

groups <- unique(groups$group)

write_parquet(sp_status,
              file.path(outf, glue::glue("metric=richness_model={acro}_what=splist.parquet")))


# Get thresholds
thresh_list <- c("p10", "mtp", "max_spec_sens")

thresholds <- purrr::map(sp_status$taxonID, function(x) {
  to_open <- file.path(
    results_folder,
    paste0("taxonid=", x, "/model=mpaeu/metrics/taxonid=", x, "_model=mpaeu_what=thresholds.parquet")
  )
  thf <- arrow::read_parquet(to_open)
  thf <- thf[, c("model", thresh_list)]
  thf$taxonID <- x
  thf
}, .progress = T)
thresholds <- bind_rows(thresholds)

get_threshold <- function(metric) {
  thresholds %>%
    select(taxonID, model, all_of(metric)) %>%
    tidyr::pivot_wider(names_from = "model", values_from = all_of(metric)) %>%
    mutate(across(2:ncol(.), ~ .x * 100))
}

thresholds_p10 <- get_threshold("p10")
thresholds_mtp <- get_threshold("mtp")
thresholds_mss <- get_threshold("max_spec_sens")


# Join models ----
models <- c("ensemble", "maxent", "rf", "xgboost", "esm", "which_valid")
scenarios <- c("current", paste0(paste0("ssp", c(126, 245, 370, 460, 585)),
                                 rep(paste0("_dec", c(50, 100)), each = 5)))
thresholds <- c("p10", "mtp", "mss")

force <- FALSE

for (i in seq_along(models)) {
  tgm <- models[i]
  for (g in seq_along(groups)) {
    tgg <- groups[g]
    for (s in seq_along(scenarios)) {
      tgs <- scenarios[s]
      for (tr in seq_along(thresholds)) {

        tgt <- thresholds[tr]

        cat("Processing model", i, "out of", length(models), "\n")
        cat("Processing group", g, "out of", length(groups), "\n")
        cat("Processing scenario", s, "out of", length(scenarios), "\n")
        cat("Processing threshold", tr, "out of", length(thresholds), "\n")
        cat("Combination - ", paste(tgm, tgg, tgs, tgt, sep = ":"), "\n\n")
        
        mtgg <- gsub('\\/', '-', tgg)
        tc_tgm <- ifelse(tgm == "which_valid", "combined", tgm)
        to_check <- tolower(file.path(outf, 
                              glue::glue(
        "metric=richness_model={acro}_method={tc_tgm}_scen={tgs}_group={mtgg}_type={tgt}_cont_cog.tif"
                              )))
        if (file.exists(to_check)) {
          if (!force) {
            cat("Already done. Skipping.\n")
            next
          }
        }
        
        sel_ids <- sp_status[sp_status$group == tgg, ]
        
        if (tgm == "which_valid") {
          sel_ids <- sel_ids[, c("taxonID", tgm)]
          tgm_ind <- sel_ids[!is.na(sel_ids$which_valid), tgm]
          sel_ids <- sel_ids[!is.na(sel_ids$which_valid), "taxonID"]
          tgm_value <- split(tgm_ind, ceiling(seq_along(sel_ids) / batch_size))
          tgm_name <- "combined"
        } else {
          sel_ids <- sel_ids[, c("taxonID", tgm)]
          sel_ids <- sel_ids[sel_ids[, tgm], "taxonID"]
          tgm_value <- tgm_name <- tgm
        }
        
        if (length(sel_ids) < 1) next
        
        batches <- split(sel_ids, ceiling(seq_along(sel_ids) / batch_size))
        
        bsizes <- unlist(lapply(batches, length))
        if (any(bsizes == 1)) {
          to_assign <- batches[[which(bsizes == 1)]]
          batches[[which(bsizes == 1)]] <- NULL
          batches[[length(batches)]] <- c(batches[[length(batches)]], to_assign)
        }
        
        if (tgt == "p10") {
          sel_th_table <- thresholds_p10
        } else if (tgt == "mtp") {
          sel_th_table <- thresholds_mtp
        } else if (tgt == "mss") {
          sel_th_table <- thresholds_mss
        } else {
          stop("Threshold not recognized.")
        }
        
        to_join <- future_map(
          seq_len(length(batches)),
          .f = proc_maps,
          batches = batches,
          model = tgm_value,
          threshold_table = sel_th_table,
          outfolder = outf,
          group = tgg,
          scenario = tgs,
          results_folder = results_folder,
          studyarea = starea,
          mask_layer = mask_target,
          save_parquet = TRUE,
          max_mem = 0.6/n_cores,
          .progress = TRUE
        )
        
        to_join <- unlist(to_join)
        
        if (is.null(to_join) || length(to_join) < 1) next
        
        all_layers <- rast(to_join)
        
        if (nlyr(all_layers) > 1) {
          sum_layers <- app(all_layers, fun = sum, na.rm = T)
        } else {
          sum_layers <- all_layers
        }
        
        sum_layers <- as.int(sum_layers)

        bin_outf <- tolower(file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm_name}_scen={tgs}_group={mtgg}_type={tgt}_bin.tif"
          )
        ))
        
        writeRaster(sum_layers, bin_outf, overwrite = T)
        
        all_layers_c <- rast(gsub("binary", "cont", to_join))
        
        if (nlyr(all_layers_c) > 1) {
          sum_layers_c <- app(all_layers_c, fun = sum, na.rm = T)
        } else {
          sum_layers_c <- all_layers_c
        }
        
        sum_layers_c <- round(sum_layers_c, 1)

        cont_outf <- tolower(file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm_name}_scen={tgs}_group={mtgg}_type={tgt}_cont.tif"
          )
        ))
        
        writeRaster(sum_layers_c, cont_outf, overwrite = T)

        obissdm::cogeo_optim(bin_outf)
        obissdm::cogeo_optim(cont_outf)

        # Aggregate parquets
        if (file.exists(gsub("\\.tif", ".parquet", to_join)[1])) {
          parquet_outf <- tolower(file.path(
            outf,
            glue::glue(
              "metric=richness_model={acro}_method={tgm_name}_scen={tgs}_group={mtgg}_type={tgt}_bin.parquet"
            )
          ))
          for (pqf in seq_along(to_join)) {
            if (pqf == 1) {
              pq_exp <- arrow::read_parquet(gsub("\\.tif", ".parquet", to_join)[pqf])
              pq_exp <- data.table::data.table(pq_exp)
            } else {
              pq_exp_pt <- arrow::read_parquet(gsub("\\.tif", ".parquet", to_join)[pqf])
              pq_exp <- merge(pq_exp, pq_exp_pt, by = c("x", "y"), all = TRUE)
              rm(pq_exp_pt)
            }
          }
          fs::file_delete(gsub("\\.tif", ".parquet", to_join))
          #ds <- arrow::open_dataset(gsub("\\.tif", ".parquet", to_join))
          pq_exp |> arrow::write_parquet(sink = parquet_outf)
          #rm(ds)
          rm(pq_exp)
        }

        fs::file_delete(to_join)
        fs::file_delete(gsub("binary", "cont", to_join))
      }
    }
  }
}

### Join full richness
cat("\n\nCreating full richness maps\n")
f <- list.files(outf, full.names = T)
f <- f[!grepl("all", f)]

for (i in seq_along(models)) {
  tgm <- models[i]
  for (s in seq_along(scenarios)) {
    tgs <- scenarios[s]
    for (tr in seq_along(thresholds)) {
      tgt <- thresholds[tr]

      cat("Processing model", i, "out of", length(models), "\n")
      cat("Processing scenario", s, "out of", length(scenarios), "\n")
      cat("Processing threshold", tr, "out of", length(thresholds), "\n")
      cat("Combination - ", paste(tgm, tgs, tgt, sep = ":"), "\n\n")

      if (tgm == "which_valid") {
        tgm <- "combined"
      }

      bin_outf <- file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm}_scen={tgs}_group=all_type={tgt}_bin.tif"
          )
        )
      
      if (file.exists(gsub("_bin", "_bin_cog", bin_outf))) {
        if (!force) {
          cat("Already done. Skipping.\n")
          next
        }
      }

      fsel <- f[grepl(tgm, f)]       # Selected model
      fsel <- fsel[grepl(tgs, fsel)] # Selected scenario
      fsel <- fsel[grepl(tgt, fsel)] # Selected threshold
      fsel <- fsel[grepl("\\.tif", fsel)] # Only rasters

      # Binary
      fsel_bin <- fsel[grepl("_bin", fsel)]
      
      r_bin <- rast(fsel_bin)
      r_bin <- sum(r_bin, na.rm = TRUE)
      r_bin <- as.int(r_bin)

      writeRaster(r_bin, bin_outf, overwrite = T)

      # Continuous
      fsel_cont <- fsel[grepl("_cont", fsel)]
      
      r_cont <- rast(fsel_cont)
      r_cont <- sum(r_cont, na.rm = TRUE)
      r_cont <- round(r_cont, 1)

      cont_outf <- file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm}_scen={tgs}_group=all_type={tgt}_cont.tif"
          )
        )

      writeRaster(r_cont, cont_outf, overwrite = T)

      obissdm::cogeo_optim(bin_outf)
      obissdm::cogeo_optim(cont_outf)
    }
  }
}



### Join full richness raw
cat("\n\nCreating full richness maps based on raw data\n")
base_rast <- "data/env/current/thetao_baseline_depthsurf_mean.tif"
study_area <- "data/shapefiles/mpa_europe_starea_v2.shp"

raw_proc <- lapply(unique(sp_status$group), raw_processing,
                   sp_list = sp_status, base_raw = base_rast, study_area = study_area,
                   output_folder = outf, acro = acro)

raw_proc <- unlist(raw_proc, use.names = F)
raw_proc <- gsub("\\.tif", "_cog.tif", raw_proc)

original_layers <- rast(raw_proc)
original_layers <- sum(original_layers, na.rm = TRUE)
original_layers <- as.int(original_layers)

aggregated_layers <- rast(gsub("original", "aggregated", raw_proc))
aggregated_layers <- sum(aggregated_layers, na.rm = TRUE)
aggregated_layers <- as.int(aggregated_layers)

outraw <- file.path(outf, glue::glue(
  "metric=richness_model={acro}_method=raw_scen=current_group=all_type=original.tif"
))
outrawagg <- file.path(outf, glue::glue(
  "metric=richness_model={acro}_method=raw_scen=current_group=all_type=aggregated.tif"
))

writeRaster(original_layers, outraw, overwrite = T)
writeRaster(aggregated_layers, outrawagg, overwrite = T)
obissdm::cogeo_optim(outraw)
obissdm::cogeo_optim(outrawagg)



### Calculate LCBD
# For later
# r <- arrow::read_parquet("https://mpaeu-dist.s3.amazonaws.com/results/diversity/metric=richness_model=mpaeu_method=combined_scen=current_group=Annelida_type=mtp_bin.parquet")

# get_betadiv <- function(div_file) {

#   con <- dbConnect(duckdb::duckdb())

#   n_count <- dbGetQuery(con, glue::glue("
#     SELECT count(*) as total
#     FROM read_parquet('{div_file}')
#   ")) |> pull()

#   batch_rows <- seq(1000, n_count, by = 1000)
#   batch_rows <- c(0, batch_rows)

#   lcbd <- numeric(n_count)

#   pb <- progress::progress_bar$new(total = length(batch_rows))

#   for (i in batch_rows) {
#     pb$tick()
#     batch <- dbGetQuery(con, glue::glue("
#       SELECT *
#       FROM read_parquet('{div_file}')
#       OFFSET {i}
#       LIMIT 1000
#     "))
#     batch <- data.table::data.table(batch)
#     batch <- data.table::setnafill(batch, fill = 0)

#     batch_coords <- batch[,1:2]
#     batch_matrix <- batch[,3:ncol(batch)]

#     beta_pred <- beta.div(batch_matrix)


#   }
  

#   ds <- arrow::open_dataset(div_file)

#   cli::cli_progress_step("Calculating Beta Diversity", spinner = T)
#   beta_div <- ds |> 
#     #slice_head(n = 500) |> 
#     select(-1, -2) |> 
#     map_batches(function(batch) {
#       batch %>%
#         as.data.frame() %>%
#         data.table::data.table() %>%
#         data.table::setnafill(fill = 0) %>%
#         mutate(lcbd = calc_metrics(.)) %>%
#         mutate(richness = rowSums(.[,1:(ncol(.)-1)])) %>%
#         select(lcbd, richness) %>%
#         as_record_batch()
#     }) |> 
#     collect()
#   cli::cli_progress_done()

#   row_total <- ds |> count() |> collect() |> pull()

#   batches <- split(seq_len(row_total), ceiling(seq_len(row_total)/1000))
#   if (length(batches[[length(batches)]]) < 2) {
#     batches[[(length(batches)-1)]] <- c(
#       batches[[(length(batches)-1)]] , batches[[length(batches)]]
#     )
#     batches[[length(batches)]] <- NULL
#   }

#   results <- lapply(seq_len(length(batches)), NULL)

#   for (bt in seq_len(length(batches))) {
#     div_matrix <- ds |> 
#       slice(batches[[bt]])
#     div_matrix <- data.table::data.table()
#     div_matrix <- data.table::setnafill(div_matrix, fill = 0)

#     div_coords <- div_matrix[,1:2]
#     div_matrix <- div_matrix[,3:ncol(div_matrix)]

#     y_sum <- rowSums(div_matrix, na.rm = T)

#     div_matrix <- div_matrix[y_sum > 0,]

#     batches <- split(seq_len(nrow(div_matrix)), ceiling(seq_len(nrow(div_matrix))/500))

#     beta_pred <- beta.div(div_matrix[1:10,])
#   }


# }



# #END
# calc_metrics <- function(div_matrix) {
#   beta.div(div_matrix)$LCBD
# }



# library(DBI)
# library(duckdb)

# # Connect to DuckDB (in-memory database for this example)
# con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

# # Query rows 100 to 1000 (using OFFSET 99 LIMIT 901)
# df_subset <- dbGetQuery(con, "
#   SELECT *
#   FROM read_parquet('path/to/file.parquet')
#   OFFSET 99
#   LIMIT 901
# ")

# df_subset