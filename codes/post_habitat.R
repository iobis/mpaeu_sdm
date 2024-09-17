############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
######################## Habitat models (Stacked SDM) ##########################

library(terra)
library(glue)
acro <- "mpaeu"

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

join_pred_sp <- function(results_folder, model_acro, species, m, sc) {

    thresholds <- c("p10", "mtp", "max_spec_sens")

    th_list <- lapply(species, function(x){
        f <- glue("{results_folder}/taxonid={x}/model={model_acro}/metrics/taxonid={x}_model={model_acro}_what=thresholds.parquet")
        arrow::read_parquet(f)
    })

    results_list <- lapply(1:length(thresholds), function(x) NULL)

    for (th in thresholds) {
        mask_file <- rast(
            file.path(
                results_folder,
                paste0(
                    glue("taxonid={species[1]}/model={model_acro}/predictions/"),
                    glue("taxonid={species[1]}_model={model_acro}_mask_cog.tif")
                )
            )
        )
        mask_file_depth <- mask_file$fit_region
        mask_file_depth[mask_file_depth == 0] <- NA
        mask_file <- mask_file$fit_ecoregions
        mask_file <- terra::mask(mask_file, mask_file_depth)
        mask_file[mask_file == 0] <- NA

        joint_raster <- rast(
            file.path(
                results_folder,
                paste0(
                    glue("taxonid={species[1]}/model={model_acro}/predictions/"),
                    glue("taxonid={species[1]}_model={model_acro}_method={m}_scen={sc}_cog.tif")
                )
            )
        )
        joint_raster <- terra::mask(joint_raster[[1]], mask_file)
        joint_raster <- joint_raster / 100

        th_val <- th_list[[1]][th_list[[1]]$model == m, th]
        joint_raster[joint_raster < as.numeric(th_val)] <- 0

        joint_raster_bin <- joint_raster
        joint_raster_bin[joint_raster_bin > 0] <- 1

        rm(mask_file, mask_file_depth)

        for (sp in 2:length(species)) {
            mask_file <- rast(
                file.path(
                    results_folder,
                    paste0(
                        glue("taxonid={species[sp]}/model={model_acro}/predictions/"),
                        glue("taxonid={species[sp]}_model={model_acro}_mask_cog.tif")
                    )
                )
            )
            mask_file_depth <- mask_file$fit_region
            mask_file_depth[mask_file_depth == 0] <- NA
            mask_file <- mask_file$fit_ecoregions
            mask_file <- terra::mask(mask_file, mask_file_depth)
            mask_file[mask_file == 0] <- NA

            species_raster <- rast(
                file.path(
                    results_folder,
                    paste0(
                        glue("taxonid={species[sp]}/model={model_acro}/predictions/"),
                        glue("taxonid={species[sp]}_model={model_acro}_method={m}_scen={sc}_cog.tif")
                    )
                )
            )
            species_raster <- terra::mask(species_raster[[1]], mask_file)
            species_raster <- species_raster / 100

            th_val <- th_list[[sp]][th_list[[sp]]$model == m, th]
            species_raster[species_raster < as.numeric(th_val)] <- 0

            species_raster_bin <- species_raster
            species_raster_bin[species_raster_bin > 0] <- 1

            joint_raster <- sum(joint_raster, species_raster, na.rm = T)
            joint_raster_bin <- sum(joint_raster_bin, species_raster_bin, na.rm = T)

            rm(mask_file, mask_file_depth, species_raster, th_val)
        }

        results_list[[which(th == thresholds)]] <- list(
            cont = joint_raster,
            bin = joint_raster_bin
        )
    }

    names(results_list) <- thresholds

    return(results_list)
}

join_predictions <- function(species, model_acro, hab_name,
                            results_folder = "results", out_folder = "habitat") {

    fs::dir_create(out_folder)

    all_preds <- list.files(file.path(results_folder,
        paste0("taxonid=", species),
        paste0("model=", model_acro), "predictions"))
    
    all_preds <- all_preds[grepl("method=", all_preds)]

    methods <- sub(".*method=([^_]+(?:_[^_]+)*?)_scen=.*", "\\1", all_preds)
    methods <- unique(methods)

    # Check exists in all
    m_ok <- lapply(methods, function(x){
        sp_cont <- lapply(species, function(sp) {
            f <- list.files(file.path(
                results_folder, paste0("taxonid=", sp, "/model=", model_acro, "/predictions/")
            ))
            any(grepl(paste0("method=", x), f))
        })
        all(unlist(sp_cont))
    })
    methods <- methods[unlist(m_ok)]

    scens <- sub(".*scen=([^_]+(?:_dec[0-9]+)?).*", "\\1", all_preds)
    scens <- unique(scens)

    cat(glue("Joining predictions for {length(species)} species"),
        glue("in {length(methods)} methods and {length(scens)} scenarios.\n"))

    for (m in methods) {
        for (sc in scens) {
            result <- join_pred_sp(results_folder, model_acro,
                species, m, sc)

            res_nams <- names(result)

            for (r in 1:length(result)) {
                outfile <- glue(
                    "habitat={hab_name}_model={model_acro}_method={m}_scen={sc}_type=cont_threshold={res_nams[r]}.tif"
                )
                outfile <- file.path(out_folder, outfile)
                trast <- result[[r]]$cont
                trast <- trast * 100
                trast <- as.int(trast)
                names(trast) <- paste0("continuous_", res_nams[r])
                writeRaster(trast, outfile, overwrite = T, datatype = "INT1U")
                obissdm::cogeo_optim(outfile)

                outfile <- glue(
                    "habitat={hab_name}_model={model_acro}_method={m}_scen={sc}_type=bin_threshold={res_nams[r]}.tif"
                )
                outfile <- file.path(out_folder, outfile)
                trast <- result[[r]]$bin
                trast <- trast * 100
                trast <- as.int(trast)
                names(trast) <- paste0("binary_", res_nams[r])
                writeRaster(trast, outfile, overwrite = T, datatype = "INT1U")
                obissdm::cogeo_optim(outfile)
            }
        }
    }
    jsonlite::write_json(
        list(habitat = hab_name,
            species = species,
            model_acro = model_acro,
            date = Sys.Date(),
            methods_used = methods,
            scenarios_used = scens),
        path = file.path(out_folder, paste0("habitat=", hab_name, "_model=", model_acro, "_what=log.json")),
        pretty = T
    )
    cat("Concluded. \n")
    return(invisible(NULL))
}

for (k in 1:length(biogenic_groups)) {
    tg_biog <- names(biogenic_groups)[k]
    cat("Processing", tg_biog, "\n")

    join_predictions(biogenic_groups[[k]], model_acro = acro, hab_name = tg_biog)
}