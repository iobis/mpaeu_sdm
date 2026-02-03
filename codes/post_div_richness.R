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
library(furrr)
source("functions/post_div_functions.R")
set.seed(2025)

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
groups <- c(groups, "all")

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
  setMinMax(r)
  mx <- max(minmax(r))
  if (mx > 65000) {
    dt <- "INT4U"
  } else {
    dt <- "INT2U"
  }
  writeRaster(r, file, datatype = dt, overwrite = TRUE)
  if (as_cog) obissdm::cogeo_optim(file)
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
      if (!is.na(decade)) {
        if (decade == 2050) {
          decade <- "dec50"
        } else if (decade == 2100) {
          decade <- "dec100"
        }
      } else {
        decade <- ""
      }
      outgroup <- gsub("\\/", "-", tolower(groups[g]))
      outfile <- glue::glue(
        "_model={model_acro}_scen={scen}{ifelse(decade == '', '', paste0('_', decade))}_",
        "th={thresh}_type={mtype}_group={outgroup}_"
      )
      if (scen == "current") {
        scen_list <- av_species |>
          filter(scenario == scen)
      } else {
        scen_list <- av_species |>
          filter(scenario == scen) |>
          filter(period == decade)
      }

      if (nrow(scen_list) != length(unique(scen_list$taxonid))) stop("Potential problem with list. Check.")

      if (interactive()) {
        layers <- rast(file.path(preproc_folder, scen_list$file))
      } else {
        layers <- vrt(file.path(preproc_folder, scen_list$file), options = c("-separate"))
      }

      # Produce traditional layer
      cont_richness <- sum(layers)

      # Produce binary layer
      bin_layers <- classify(layers, matrix(c(0, Inf, 1), nrow = 1))
      bin_richness <- sum(bin_layers)

      # SESAM probabilities
      sesam_richness <- sesam_prr(layers, (cont_richness/100))
      sesam_richness <- sum(sesam_richness)

      # Save rasters
      base_f <- paste0("metric=richness", outfile)
      tempf <- paste0(base_f, "what=continuous.tif")
      save_raster(cont_richness, file.path(out_folder, tempf))
      tempf <- paste0(base_f, "what=binary.tif")
      save_raster(bin_richness, file.path(out_folder, tempf))
      tempf <- paste0(base_f, "what=sesam.tif")
      save_raster(sesam_richness, file.path(out_folder, tempf))

      # LCBD
      coarse_binary <- aggregate(bin_layers, 5, max, na.rm = T)
      coarse_binary_df <- as.data.frame(coarse_binary, cell = T)
      cells <- coarse_binary_df[, 1]
      coarse_binary_df <- coarse_binary_df[, 2:ncol(coarse_binary_df)]
      valid_cells <- rowSums(coarse_binary_df)
      valid_cells <- which(valid_cells > 0)
      cells <- cells[valid_cells]
      coarse_binary_df <- coarse_binary_df[valid_cells, ]
      beta_div <- lcbd_cpp(as.matrix(coarse_binary_df))
      lcbd_rast <- coarse_binary[[1]]
      lcbd_rast[!is.na(lcbd_rast)] <- 0
      beta_div <- beta_div * (1 / max(beta_div))
      lcbd_rast[cells] <- beta_div
      lcbd_rast <- lcbd_rast * 100

      base_f <- paste0("metric=lcbd", outfile)
      tempf <- paste0(gsub("_$", "", base_f), ".tif")
      save_raster(lcbd_rast, file.path(out_folder, tempf))
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

    tempf <- paste0("metric=richness_model=", model_acro, "_group=", gsub("\\/", "-", tolower(groups[g])), "_what=raw.tif")
    save_raster(true_richness, file.path(out_folder, tempf))
  }
}



# Step 3 - save species list -----
colnames(proc_list)[1] <- "taxonID"
proc_list$taxonID <- as.integer(proc_list$taxonID)
proc_list <- left_join(proc_list, species_list[,c("taxonID", "scientificName", "group")])

# Check validity on study area
pb <- progress::progress_bar$new(total = nrow(proc_list))
proc_list$valid_on_area <- FALSE
for (k in seq_len(nrow(proc_list))) {
  pb$tick()
  r <- rast(file.path(preproc_folder, proc_list$file[k]))
  if (max(minmax(r)) > 0) {
    proc_list$valid_on_area[k] <- TRUE
  }
}

proc_list_dist <- proc_list |>
  filter(valid_on_area) |>
  select(-period, -scenario, -file, -valid_on_area) |>
  group_by(taxonID, scientificName, group, method, threshold, type) |>
  distinct()

arrow::write_parquet(proc_list_dist, file.path(out_folder, glue::glue("metric=richness_model={model_acro}_what=splist.parquet")))



# Step 4 - refugia analysis -------
parallel <- TRUE
if (parallel) {
  n_workers <- 10
  plan(multisession, workers = n_workers)
  max_mem <- 0.9 / n_workers
}

blank_raster <- tempfile(fileext = ".tif")
blank_m <- rast(file.path(preproc_folder, proc_list$file[1]))
blank_m[!is.na(blank_m)] <- 0
writeRaster(blank_m, blank_raster)

name_and_save <- \(r, outfile, what) {
  base_f <- paste0("metric=refugia", outfile)
  tempf <- paste0(base_f, "what=", what, ".tif")
  save_raster(as.int(r), file.path(out_folder, tempf))
  return(invisible())
} 

message("=========== Starting refugia analysis ===========\n\n")
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

    # If not running in parallel, use below:
    # Prepare current layer -----
    current_list <- av_species |>
            filter(scenario == "current")
    if (!parallel) {
      
      if (interactive()) {
        current_rast <- rast(file.path(preproc_folder, current_list$file))
      } else {
        current_rast <- vrt(file.path(preproc_folder, current_list$file), options = c("-separate"))
      }

      current_rast <- classify(current_rast, matrix(c(0, Inf, 1), ncol = 3))

      richness <- sum(current_rast)
    }

    # Process future layers -----
    future_map(2:nrow(scen_grid), \(sc, model_acro = model_acro) {
      require(terra)
      terraOptions(memfrac = max_mem)
      scen <- scen_grid$scenario[sc]
      decade <- scen_grid$decade[sc]
      if (decade == 2050) {
        decade <- "dec50"
      } else if (decade == 2100) {
        decade <- "dec100"
      }
      outgroup <- gsub("\\/", "-", tolower(groups[g]))
      outfile <- glue::glue(
        "_model={model_acro}_scen={scen}{ifelse(decade == '', '', paste0('_', decade))}_",
        "th={thresh}_type={mtype}_group={outgroup}_"
      )
      scen_list <- av_species |>
          filter(scenario == scen) |>
          filter(period == decade)

      if (nrow(scen_list) != length(unique(scen_list$taxonid))) stop("Potential problem with list. Check.")

      scen_list <- scen_list |>
        select(taxonid, future_file = file) |>
        right_join(current_list) |>
        mutate(future_file = ifelse(is.na(future_file), blank_raster, future_file)) |>
        select(taxonid, file = future_file)
      
      if (nrow(scen_list) != nrow(current_list)) stop("Potential problem with list. Check.")

      if (parallel) {
        
        if (interactive()) {
          current_rast <- rast(file.path(preproc_folder, current_list$file))
        } else {
          current_rast <- vrt(file.path(preproc_folder, current_list$file), options = c("-separate"))
        }

        current_rast <- classify(current_rast, matrix(c(0, Inf, 1), ncol = 3))

        richness <- sum(current_rast)
      }

      if (interactive()) {
        layers <- rast(file.path(preproc_folder, scen_list$file))
      } else {
        layers <- vrt(file.path(preproc_folder, scen_list$file), options = c("-separate"))
      }

      layers <- classify(layers, matrix(c(0, Inf, 3), ncol = 3))
      # 3 - 1 ==> 2 = stability = keep present
      # 3 - 0 ==> 3 = gain
      # 0 - 1 ==> -1 = loss
      # 0 - 0 ==> 0 = stability = keep absent

      result <- layers - current_rast

      gain <- app(result, \(x) sum(x == 3))
      loss <- app(result, \(x) sum(x == -1))

      stability_pres <- app(result, \(x) sum(x == 2))
      stability_absc <- app(result, \(x) sum(x == 0))

      turnover <- (loss + gain)/(richness + gain)

      refugia <- 100 - ((loss * 100)/richness)

      gain |> name_and_save(outfile, "gain")
      loss |> name_and_save(outfile, "loss")
      stability_pres |> name_and_save(outfile, "stability_pres")
      stability_absc |> name_and_save(outfile, "stability_absc")
      (turnover * 100) |> name_and_save(outfile, "turnover")
      (refugia * 100) |> name_and_save(outfile, "refugia")

      return(invisible(NULL))
    }, .progress = TRUE)
  }
}