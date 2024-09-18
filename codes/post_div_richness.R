############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
#################### Diversity maps - richness and others ######################

# Load packages/settings -----
library(terra)
library(glue)
library(future)
library(furrr)
acro <- "mpaeu"
results_folder <- "results"
output_folder <- "diversity"
thresh_list <- c("p10", "mtp", "max_spec_sens")
fs::dir_create(output_folder)


# List models to assembly ----
done_sp <- list.files(results_folder)

done_good <- lapply(done_sp, function(x) {
    f <- list.files(file.path(results_folder, x, paste0("/model=", acro)))
    any(grepl("log.json", f))
})

done_good <- done_sp[unlist(done_good)]
done_good <- as.numeric(gsub("taxonid=", "", done_good))

methods_within <- lapply(done_good, function(sp) {

    logf <- jsonlite::read_json(
        file.path(
            results_folder,
            paste0("taxonid=", sp),
            paste0("model=", acro), paste0("taxonid=", sp, "_model=", acro, "_what=log.json")
    ))

    success <- unlist(logf$model_good)

    if (success[[1]] == 1) {
        success <- "esm"
    }

    if (length(success) > 1 && !is.null(logf$model_result$ensemble[[1]])) {
        success <- c(success, "ensemble")
    }
    success
})

methods_within <- as.data.frame(table(unlist(methods_within)))
max_val <- max(methods_within$Freq)
methods_within <- methods_within[methods_within$Freq == max_val,]

methods <- as.character(methods_within[,1])

all_preds <- list.files(file.path(results_folder, paste0("taxonid=", done_good[1]),
    paste0("model=", acro), "predictions"))
all_preds <- all_preds[!grepl("mask|mess|shape", all_preds)]
scens <- sub(".*scen=([^_]+(?:_dec[0-9]+)?).*", "\\1", all_preds)
scens <- unique(scens)

cli::cli_alert_info("Joining predictions for {length(done_good)} species in {length(methods)} models and {length(scens)} scenarios.")


mask_thresh <- function(species, model_acro, m, sc, th_val = NULL, mask_layer) {

    species_raster <- rast(
        file.path(
            results_folder,
            paste0(
                glue("taxonid={species[1]}/model={model_acro}/predictions/"),
                glue("taxonid={species[1]}_model={model_acro}_method={m}_scen={sc}_cog.tif")
            )
        )
    )
    species_raster <- terra::mask(species_raster[[1]], mask_layer)
    species_raster <- species_raster / 100

    if (!is.null(th_val)) {
        species_raster[species_raster < as.numeric(th_val)] <- 0

        species_raster_bin <- species_raster
        species_raster_bin[species_raster_bin > 0] <- 1

        result <- list(cont = species_raster, bin = species_raster_bin)
    } else {
        result <- species_raster
    }

    return(result)
}

# TODO: add methods from Zurrell/Calabrese to adjust predictions
join_models <- function(species, thresholds = c("p10", "mtp", "max_spec_sens"),
                        model_acro, m, sc, results_folder) {
    modes <- c("continuous", paste0("thresh_", thresholds), paste0("bin_", thresholds))

    results <- lapply(1:length(modes), function(x) NULL)
    names(results) <- modes

    for (i in 1:length(species)) {
        example_file <- file.path(results_folder,
                        glue("taxonid={species[i]}/model={model_acro}/predictions/"),
                        glue("taxonid={species[i]}_model={model_acro}_method={m}_scen={sc}_cog.tif")
                    )

        if (file.exists(example_file)) {
            mask_file <- rast(
                file.path(
                    results_folder,
                    paste0(
                        glue("taxonid={species[i]}/model={model_acro}/predictions/"),
                        glue("taxonid={species[i]}_model={model_acro}_mask_cog.tif")
                    )
                )
            )
            mask_file_depth <- mask_file$fit_region
            mask_file_depth[mask_file_depth == 0] <- NA
            mask_file <- mask_file$fit_ecoregions
            mask_file <- terra::mask(mask_file, mask_file_depth)
            mask_file[mask_file == 0] <- NA

            if (is.null(results$continuous)) {
                results$continuous <- mask_thresh(
                    species = species[i], model_acro = model_acro, m = m, sc = sc,
                    th_val = NULL, mask_layer = mask_file
                )
            } else {
                results$continuous <- sum(results$continuous, mask_thresh(
                    species = species[i], model_acro = model_acro, m = m, sc = sc,
                    th_val = NULL, mask_layer = mask_file
                ), na.rm = T)
            }

            th_file <- arrow::read_parquet(file.path(
                glue("{results_folder}/taxonid={species[i]}/model={model_acro}/"),
                glue("metrics/taxonid={species[i]}_model={model_acro}_what=thresholds.parquet")
            ))

            for (th in 1:length(thresholds)) {
                th_val <- as.numeric(th_file[th_file$model == m, thresholds[th]])

                ready_rast <- mask_thresh(
                    species = species[i], model_acro = model_acro, m = m, sc = sc,
                    th_val = th_val, mask_layer = mask_file
                )

                thc_n <- paste0("thresh_", thresholds[th])
                thb_n <- paste0("bin_", thresholds[th])

                if (is.null(results[[thc_n]])) {
                    results[[thc_n]] <- ready_rast$cont
                    results[[thb_n]] <- ready_rast$bin
                } else {
                    results[[thc_n]] <- sum(results[[thc_n]], ready_rast$cont, na.rm = T)
                    results[[thb_n]] <- sum(results[[thb_n]], ready_rast$bin, na.rm = T)
                }
            }
        }
    }

    return(results)
}

# Get richness ----
to_get <- expand.grid(methods = methods, scens = scens)

plan(multisession, workers = 5)

future_map(1:nrow(to_get), function(id){
    tm <- to_get$methods[id]
    ts <- to_get$scens[id]

    results_list <- join_models(
        species = done_good,
        thresholds = thresh_list,
        model_acro = acro, m = tm, sc = ts, results_folder = results_folder
    )

    results_names <- names(results_list)
    results_names <- sub("_", "_threshold=", results_names)

    r <- lapply(1:length(results_list), function(x) {
        outfile <- glue(
            "metric=richness_model={acro}_method={tm}_scen={ts}_type={results_names[x]}.tif"
        )
        writeRaster(results_list[[x]], file.path(output_folder, outfile), overwrite = T)
        obissdm::cogeo_optim(file.path(output_folder, outfile))
    })

    return(invisible(NULL))
}, .progress = T)

jsonlite::write_json(
        list(metrics = "richness",
            species = done_good,
            model_acro = model_acro,
            date = Sys.Date(),
            methods_used = methods,
            scenarios_used = scens,
            thresholds = thresh_list),
        path = file.path(output_folder, paste0("metric=richness", "_model=", model_acro, "_what=log.json")),
        pretty = T
)

cli::cli_alert_success("Files saved in {.path {output_folder}}")
