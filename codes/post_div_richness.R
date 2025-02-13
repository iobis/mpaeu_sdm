############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
#################### Diversity maps - richness and others ######################

# Load packages/settings -----
library(progress)
library(arrow)
library(dplyr)
library(terra)
library(furrr)
source("functions/post_div_functions.R")



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
        to_check <- file.path(outf, 
                              glue::glue(
        "metric=richness_model={acro}_method={tgm}_scen={tgs}_group={mtgg}_type={tgt}_cont_cog.tif"
                              ))
        if (file.exists(to_check)) {
          if (!force) {
            cat("Already done. Skipping.\n")
            next
          }
        }
        
        sel_ids <- sp_status[sp_status$group == tgg, ]
        
        if (tgm == "which_valid") {
          sel_ids <- sel_ids[, c("taxonID", tgm)]
          tgm <- sel_ids[!is.na(sel_ids$which_valid), tgm]
          sel_ids <- sel_ids[!is.na(sel_ids$which_valid), "taxonID"]
          tgm <- split(tgm, ceiling(seq_along(sel_ids) / batch_size))
        } else {
          sel_ids <- sel_ids[, c("taxonID", tgm)]
          sel_ids <- sel_ids[sel_ids[, tgm], "taxonID"]
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
          model = tgm,
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

        bin_outf <- file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm}_scen={tgs}_group={mtgg}_type={tgt}_bin.tif"
          )
        )
        
        writeRaster(sum_layers, bin_outf, overwrite = T)
        
        all_layers_c <- rast(gsub("binary", "cont", to_join))
        
        if (nlyr(all_layers_c) > 1) {
          sum_layers_c <- app(all_layers_c, fun = sum, na.rm = T)
        } else {
          sum_layers_c <- all_layers_c
        }
        
        sum_layers_c <- round(sum_layers_c, 1)

        cont_outf <- file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm}_scen={tgs}_group={mtgg}_type={tgt}_cont.tif"
          )
        )
        
        writeRaster(sum_layers_c, cont_outf, overwrite = T)

        obissdm::cogeo_optim(bin_outf)
        obissdm::cogeo_optim(cont_outf)

        # Aggregate parquets
        if (file.exists(gsub("\\.tif", ".parquet", to_join)[1])) {
          parquet_outf <- file.path(
            outf,
            glue::glue(
              "metric=richness_model={acro}_method={tgm}_scen={tgs}_group={mtgg}_type={tgt}_bin.parquet"
            )
          )
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

      fsel <- f[grepl(tgm, f)]       # Selected model
      fsel <- fsel[grepl(tgs, fsel)] # Selected scenario
      fsel <- fsel[grepl(tgt, fsel)] # Selected threshold

      # Binary
      fsel_bin <- fsel[grepl("_bin", fsel)]
      
      r_bin <- rast(fsel_bin)
      r_bin <- sum(r_bin)
      r_bin <- as.int(r_bin)

      bin_outf <- file.path(
          outf,
          glue::glue(
            "metric=richness_model={acro}_method={tgm}_scen={tgs}_group=all_type={tgt}_bin.tif"
          )
        )

      writeRaster(r_bin, bin_outf, overwrite = T)

      # Continuous
      fsel_cont <- fsel[grepl("_cont", fsel)]
      
      r_cont <- rast(fsel_cont)
      r_cont <- sum(r_cont)
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

#END
