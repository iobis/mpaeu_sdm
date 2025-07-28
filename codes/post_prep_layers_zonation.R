############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# July of 2025
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################## Post-processing #################################
####################### Prepare layers for Zonation ############################

library(glue)
library(dplyr)
library(terra)

# Settings --------
out_folder <- "/data/scps/processed_layers_v2"
fs::dir_create(out_folder)
results_folder <- "/data/scps/v5/results"
parallel <- FALSE
n_cores <- 80
max_mem <- (0.9/n_cores)

# List available species
species <- list.files(results_folder)
species <- gsub("taxonid=", "", species)

# Create a function to check good model -------
check_model <- function(results_folder, sp) {
    l <- jsonlite::read_json(glue(
        "{results_folder}/taxonid={sp}/model=mpaeu/taxonid={sp}_model=mpaeu_what=log.json"
    ))
    goodm <- unlist(l$model_good, use.names = F)
    if (is.numeric(goodm)) {
        goodm <- "esm"
    }
    if ("ensemble" %in% goodm) {
        best_m <- "ensemble"
    } else {
        evals <- lapply(goodm, \(mm) {
            mmb <- ifelse(mm == "rf", "rf_classification_ds", mm)
            eval <- arrow::read_parquet(glue(
                "{results_folder}/taxonid={sp}/model=mpaeu/metrics/taxonid={sp}_model=mpaeu_method={mmb}_what=cvmetrics.parquet"
            ))
            eval <- eval |>
                select(AUC = auc, CBI = cbi, TSS = tss_maxsss)
            eval <- round(apply(eval, 2, mean, na.rm = T), 2)
            eval <- as.data.frame(t(eval))
            eval$model <- mm
            eval
        })
        evals <- bind_rows(evals)

        best_m <- evals$model[which.max(evals$CBI)]
        best_m <- ifelse(best_m == "rf", "rf_classification_ds", best_m)
    }
    return(best_m)
}

# Create a function to do the processing --------
proc_layers <- function(sp, results_folder, out_folder, global_mask,
                        sel_threshold = "p10", type = "std",
                        alpha = 0.5, model_acro = "mpaeu",
                        check_exists = TRUE, verbose = TRUE, max_mem = NULL) {
    if (verbose) message("Processing ", sp, " #", which(species == sp))

    best_m <- check_model(results_folder, sp)

    if (check_exists) {
        if (file.exists(file.path(
            out_folder,
            glue("taxonid={sp}_model={model_acro}_method={best_m}_scen=current_th={sel_threshold}_type={type}.tif")
        ))) {
            return(paste0(sp, "_done"))
        }
    }

    global_mask <- vect(global_mask)

    # Load mask
    masks <- rast(glue(
        "{results_folder}/taxonid={sp}/model={model_acro}/predictions/taxonid={sp}_model={model_acro}_mask_cog.tif"
    ))
    masks <- masks$fit_region
    NAflag(masks) <- 0

    # Load threshold
    thresholds <- arrow::read_parquet(glue(
        "{results_folder}/taxonid={sp}/model={model_acro}/metrics/taxonid={sp}_model={model_acro}_what=thresholds.parquet"
    ))

    th <- thresholds[grepl(strtrim(best_m, 2), thresholds$model), ]
    if (sel_threshold == "mss") {
        th <- th |> pull(max_spec_sens) * 100
    } else if (sel_threshold == "p10") {
        th <- th |> pull(p10) * 100
    } else {
        stop("Threshold not available")
    }

    all_preds <- list.files(glue("{results_folder}/taxonid={sp}/model={model_acro}/predictions"), full.names = T)
    all_preds <- all_preds[grepl(best_m, all_preds)]

    model_preds <- all_preds[!grepl("what=bootcv", all_preds)]
    boot_preds <- all_preds[grepl("what=bootcv", all_preds)]

    if (length(model_preds) != length(boot_preds)) {
        return(paste0(sp, "_boot-nav-for-all"))
    } else if (length(boot_preds) == 0) {
        return(paste0(sp, "_no-boot-files"))
    }

    base <- model_preds[1] |>
        rast() |>
        terra::classify(matrix(data = c(-Inf, Inf, 0), nrow = 1), right = FALSE) |>
        crop(global_mask) |>
        mask(global_mask)

    for (lp in model_preds) {
        lyr_pred <- rast(lp)[[1]]
        lyr_boot <- rast(gsub("_cog.tif", "_what=bootcv_cog.tif", lp))[["sd"]]

        lyr_boot <- lyr_boot * alpha

        lyr_pred <- terra::classify(lyr_pred, matrix(data = c(-Inf, th, 0), nrow = 1), right = FALSE)

        lyr_pred <- lyr_pred - lyr_boot

        lyr_pred <- lyr_pred |>
            mask(masks) |>
            mask(global_mask) |>
            crop(global_mask) |>
            sum(base, na.rm = TRUE)

        if (minmax(lyr_pred)["min",] < 0) {
            lyr_pred <- terra::classify(lyr_pred, matrix(data = c(-Inf, 0, 0), nrow = 1), right = FALSE)
        }

        basef <- basename(lp)
        basef <- gsub("_cog.tif", glue("_th={sel_threshold}_type={type}.tif"), basef)
        writeRaster(
            lyr_pred,
            file.path(out_folder, basef),
            datatype = "INT1U",
            overwrite = TRUE
        )
    }

    return(paste0(sp, "_done"))
}

# Apply processing -----
if (parallel) {
    require("furrr")
    plan(multisession, workers = n_cores)

    proc_result <- future_map(species, proc_layers,
        results_folder, out_folder, global_mask,
        sel_threshold = "p10", type = "std",
        alpha = 0.5, model_acro = "mpaeu",
        check_exists = TRUE, verbose = FALSE,
        max_mem = max_mem, .progress = TRUE
    )

} else {
    proc_result <- lapply(species, proc_layers,
        results_folder, out_folder, global_mask,
        sel_threshold = "p10", type = "std",
        alpha = 0.5, model_acro = "mpaeu",
        check_exists = TRUE, verbose = TRUE
    )
}


# END
