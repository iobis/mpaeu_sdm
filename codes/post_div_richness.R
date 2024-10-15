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
library(dplyr)
library(foreach)
source("functions/post_div_functions.R")

acro <- "mpaeu"
results_folder <- "/data/scps/v3/results/"
output_folder <- "diversity"
thresh_list <- c("p10", "mtp", "max_spec_sens")
fs::dir_create(output_folder)
masks_temp <- "masks_temp"
fs::dir_create(masks_temp)

# Load shapefiles
eez <- vect("data/shapefiles/eez_v12.gpkg")
wdpa <- vect("data/shapefiles/EMODnet_HA_Environment_WDPA_ Sep2023_20230911.gdb/")

# Prepare groups of interest
sp_list <- get_sp_list()

# List models to assembly ----
done_sp <- list.files(results_folder)

done_good <- lapply(done_sp, function(x) {
    f <- list.files(file.path(results_folder, x, paste0("/model=", acro)))
    any(grepl("log.json", f))
})

done_good <- done_sp[unlist(done_good)]
done_good <- as.numeric(gsub("taxonid=", "", done_good))

# Get settings ----
all_preds <- list.files(file.path(
    results_folder, paste0("taxonid=", done_good[1]),
    paste0("model=", acro), "predictions"
))
all_preds <- all_preds[!grepl("mask|mess|shape", all_preds)]
scens <- sub(".*scen=([^_]+(?:_dec[0-9]+)?).*", "\\1", all_preds)
scens <- unique(scens)

groups <- unique(sp_list$group)

methods <- c("ensemble", "maxent", "rf", "xgboost", "raw")

groups_methods <- lapply(seq_along(methods), function(x) NULL)
names(groups_methods) <- methods
groups_sp <- lapply(seq_along(groups), function(x) groups_methods)
names(groups_sp) <- groups

base_raw <- terra::rast("data/env/current/thetao_baseline_depthsurf_mean.tif")
base_raw[!is.na(base_raw)] <- 0

cli::cli_alert_info("Joining predictions for {length(done_good)} species in {length(methods)} models and {length(scens)} scenarios.")

for (gr in groups) {
    cli::cli_alert_info("Joining predictions for group {gr}")

    sel_species <- done_good[done_good %in% sp_list$taxonID[sp_list$group == gr]]

    if (length(sel_species) < 1) {
        cli::cli_alert_danger("No species available for group {gr}. Skipping.")
        next
    }

    cli::cli_alert_info("{length(sel_species)} species available for model based join.")

    raw_processing(sp_list, gr, base_raw, output_folder, eez, wdpa)

    cli::cli_alert_info("Done. Processing modeled...")

    # Get metrics
    metrics <- lapply(sel_species, function(x) {
        m <- arrow::read_parquet(paste0(
            results_folder, "/taxonid=", x, "/model=mpaeu/metrics/taxonid=",
            x, "_model=mpaeu_what=thresholds.parquet"
        ))
        m[, c("model", thresh_list)]
    })

    # Load masks
    mask_names <- paste0(
        results_folder, "/taxonid=", sel_species, "/model=mpaeu/predictions/taxonid=",
        sel_species, "_model=mpaeu_mask_cog.tif"
    )

    cl <- parallel::makeCluster(30)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, "mask_names")
    r <- foreach(i = seq_along(mask_names)) %dopar% {
        mask_tnam <- file.path(masks_temp, basename(mask_names[i]))
        if (!file.exists(mask_tnam)) {
            mask_native <- terra::rast(mask_names[i], lyrs = 1)
            mask_fit <- terra::rast(mask_names[i], lyrs = 3)

            terra::NAflag(mask_native) <- 0
            terra::NAflag(mask_fit) <- 0

            terra::mask(mask_native, mask_fit, filename = mask_tnam, overwrite = T)
        }
        "done"
    }

    new_masks <- file.path(masks_temp, basename(mask_names[1]))

    for (me in methods) {
        cli::cli_alert_info("Joining method {me}.")

        me_ed <- ifelse(me == "rf", "rf_classification_ds", me)

        # Check which/how many species available for this method
        me_flist <- paste0(
            results_folder, "/taxonid=", sel_species, "/",
            glue("model={acro}/predictions/taxonid={sel_species}_model=mpaeu_method={me_ed}_scen=current_cog.tif")
        )
        me_flist <- me_flist[file.exists(me_flist)]

        if (length(me_flist) < 1) {
            groups_sp[[gr]][[me]] <- NA
            next
        } else {
            groups_sp[[gr]][[me]] <- basename(me_flist)
            groups_sp[[gr]][[me]] <- gsub(".*(taxonid=[0-9]+).*", "\\1", groups_sp[[gr]][[me]])
            groups_sp[[gr]][[me]] <- as.numeric(gsub("taxonid=", "", groups_sp[[gr]][[me]]))
        }

        p10 <- unlist(lapply(metrics, function(x) x$p10[x$model == me] * 100))
        mtp <- unlist(lapply(metrics, function(x) x$mtp[x$model == me] * 100))
        mss <- unlist(lapply(metrics, function(x) x$max_spec_sens[x$model == me] * 100))

        for (sc in scens) {
            cli::cli_alert_info("Joining scenario {sc}.")

            me_pred_sc <- paste0(
                results_folder, "/taxonid=", sel_species, "/",
                glue("model={acro}/predictions/taxonid={sel_species}_model=mpaeu_method={me_ed}_scen={sc}_cog.tif")
            )

            parallel::clusterExport(cl, "new_masks")
            parallel::clusterExport(cl, "me_pred_sc")
            parallel::clusterExport(cl, "thresh_list")

            r <- foreach(i = seq_along(me_pred_sc)) %dopar% {
                if (file.exists(me_pred_sc[i])) {
                    species_rast <- terra::rast(me_pred_sc[i])
                    taxonid <- gsub(".*(taxonid=[0-9]+).*", "\\1", me_pred_sc[i])
                    mnam <- gsub("taxonid=[0-9]+", taxonid, new_masks)
                    species_mask <- terra::rast(mnam)

                    sp_masked <- terra::mask(species_rast, species_mask)

                    basef <- basename(me_pred_sc[i])
                    based <- dirname(mnam)

                    for (th in thresh_list) {
                        temp_rast <- sp_masked
                        if (th == "p10") {
                            metric <- p10[i]
                        } else if (th == "mtp") {
                            metric <- mtp[i]
                        } else if (th == "max_spec_sens") {
                            metric <- mss[i]
                        }
                        temp_rast[temp_rast < metric] <- 0
                        outfile <- file.path(based, gsub("cog", paste0("type=", th, "_ct"), basef))
                        terra::writeRaster(temp_rast, outfile, overwrite = T)
                        temp_rast[temp_rast > 0] <- 1
                        outfile <- file.path(based, gsub("cog", paste0("type=", th, "_bin"), basef))
                        terra::writeRaster(temp_rast, outfile, overwrite = T)
                    }
                }

                "done"
            }

            for (th in thresh_list) {
                basef <- basename(me_pred_sc)
                cont <- gsub("cog", paste0("type=", th, "_ct"), basef)
                bin <- gsub("cog", paste0("type=", th, "_bin"), basef)
                based <- dirname(new_masks)
                cont <- file.path(based, cont)
                bin <- file.path(based, bin)

                cont <- cont[file.exists(cont)]
                bin <- bin[file.exists(bin)]

                if (length(cont) > 0) {
                    outcont <- glue(
                        "metric=richness_model={acro}_method={me}_scen={sc}_group={gsub('/', '-', tolower(gr))}_type={th}_cont.tif"
                    )
                    outbin <- glue(
                        "metric=richness_model={acro}_method={me}_scen={sc}_group={gsub('/', '-', tolower(gr))}_type={th}_bin.tif"
                    )

                    cont_rast <- rast(cont)
                    cont_rast <- sum(cont_rast, na.rm = T)

                    bin_rast <- rast(bin)
                    bin_rast <- sum(bin_rast, na.rm = T)

                    cont_rast <- as.int(cont_rast)
                    bin_rast <- as.int(bin_rast)

                    if (terra::minmax(cont_rast)[2, 1] <= 255) {
                        format_out <- "INT1U"
                    } else if (terra::minmax(cont_rast)[2, 1] <= 65535) {
                        format_out <- "INT2U"
                    } else {
                        format_out <- "INT4U"
                    }

                    terra::writeRaster(cont_rast, file.path(output_folder, outcont), overwrite = T, datatype = format_out)
                    obissdm::cogeo_optim(file.path(output_folder, outcont))

                    if (terra::minmax(bin_rast)[2, 1] <= 255) {
                        format_out <- "INT1U"
                    } else if (terra::minmax(bin_rast)[2, 1] <= 65535) {
                        format_out <- "INT2U"
                    } else {
                        format_out <- "INT4U"
                    }

                    terra::writeRaster(bin_rast, file.path(output_folder, outbin), overwrite = T, datatype = format_out)
                    obissdm::cogeo_optim(file.path(output_folder, outbin))
                }
            }
        }
    }
    cli::cli_alert_success("Group {cli::bg_cyan(gr)} concluded, removing temporary files.")
    fs::file_delete(list.files(masks_temp, full.names = T))
}

jsonlite::write_json(
    list(
        metrics = "richness",
        species = groups_sp,
        model_acro = acro,
        date = Sys.Date(),
        methods_used = methods,
        scenarios_used = scens,
        groups_used = groups,
        thresholds = thresh_list
    ),
    path = file.path(output_folder, paste0("metric=richness", "_model=", acro, "_what=log.json")),
    pretty = T
)

cli::cli_alert_success("Files saved in {.path {output_folder}}")
