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
