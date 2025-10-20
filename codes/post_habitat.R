############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
######################## Habitat models (Stacked SDM) ##########################

library(terra)
library(glue)
source("functions/post_div_functions.R")

# Settings
results_folder <- "/data/scps/v5/results"
out_folder <- "/data/scps/v5/habitat"
preproc_folder <- file.path(out_folder, "preproc")
fs::dir_create(out_folder)
fs::dir_create(preproc_folder)
model_acro <- "mpaeu"
target_thresholds <- c("p10", "mss")
target_types <- c("std", "const")

biogenic_groups <- list(
    seagrass = c(
        145794, 145793, 145795
    ),
    kelp = c(
        145725, 145728, 234483, 145735
    ),
    corals = c(
        1245747, 135209, 135168
    ),
    bivalves_beds = c(
        140480, 140481, 140467, 140658
    ),
    maerl = c(
        145165, 145170, 145199
    ),
    polychaete_reefs = c(
        130866, 130867
    )
)

# Step 1 - pre-process layers
prepare_layers(
    thresholds = target_thresholds,
    type = target_types,
    results_folder = results_folder,
    out_folder = preproc_folder,
    parallel = TRUE,
    n_cores = length(unlist(biogenic_groups, use.names = F)),
    max_mem = TRUE,
    global_mask = "data/shapefiles/mpa_europe_starea_v3.gpkg",
    base_file = "data/env/current/thetao_baseline_depthsurf_mean.tif",
    species = unlist(biogenic_groups, use.names = FALSE),
    verbose = TRUE
)

# Step 2 - assemble
types_grid <- expand.grid(
    sel_threshold = target_thresholds,
    type = target_types
)

scen_grid <- data.frame(
    scenario = c("current", rep(paste0("ssp", c(126, 245, 370, 460, 585)), 2)),
    decade = c(NA, rep(c(2050, 2100), each = 5))
)

groups <- names(biogenic_groups)

proc_list <- list_processed(preproc_folder, write_csv = FALSE)

save_raster <- \(r, file, as_cog = TRUE) {
  r <- as.int(r)
  writeRaster(r, file, datatype = "INT2U", overwrite = TRUE)
  obissdm::cogeo_optim(file)
  return(invisible(NULL))
}

for (tg in seq_len(nrow(types_grid))) {
    cli::cli_alert_info(cli::bg_cyan("Combination {tg} out of {nrow(types_grid)}"))
    thresh <- types_grid$sel_threshold[tg]
    mtype <- types_grid$type[tg]

    for (g in seq_len(length(groups))) {
        cli::cli_alert_info(cli::bg_green("Group {g} out of {length(groups)}"))

        sel_species <- biogenic_groups[[g]]

        av_species <- proc_list |>
            filter(taxonid %in% sel_species) |>
            filter(threshold == thresh) |>
            filter(type == mtype)

        if (length(unique(av_species$taxonid)) != length(sel_species)) stop("Problem with group ", groups[g])

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
                "th={thresh}_type={mtype}_"
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

            # Save rasters
            base_f <- paste0("habitat=", outgroup, outfile)
            tempf <- paste0(base_f, "what=continuous.tif")
            save_raster(cont_richness, file.path(out_folder, tempf))
            tempf <- paste0(base_f, "what=binary.tif")
            save_raster(bin_richness, file.path(out_folder, tempf))
        }
    }
}

# Save occurrence records info ------
species_list <- read.csv(obissdm::recent_file("data", "all_splist"))
species_list <- species_list[species_list$taxonID %in% unlist(biogenic_groups, use.names = F),]

for (g in seq_len(length(groups))) {
    cli::cli_alert_info(cli::bg_green("Group {g} out of {length(groups)}"))

    sel_species <- biogenic_groups[[g]]

    species_recs <- file.path(results_folder, glue(
        "taxonid={sel_species}/model={model_acro}/taxonid={sel_species}_model={model_acro}_what=fitocc.parquet"
    ))

    group_points <- vector("list", length(species_recs))
    for (ss in seq_along(sel_species)) {
        group_points[[ss]] <- arrow::read_parquet(species_recs[ss])
        group_points[[ss]]$species <- species_list$scientificName[species_list$taxonID == sel_species[ss]]
        group_points[[ss]]$taxonID <- sel_species[ss]
    }
    group_points <- dplyr::bind_rows(group_points)

    outgroup <- gsub("\\/", "-", tolower(groups[g]))
    outfile <- glue::glue(
        "_model={model_acro}_what=fitocc.parquet"
    )
    arrow::write_parquet(
        group_points, file.path(out_folder, paste0("habitat=", outgroup, outfile))
    )
}