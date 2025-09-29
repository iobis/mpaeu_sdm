library(reticulate)
library(glue)
source("functions/utils.R")
#py_install("pystac")
pystac <- import("pystac")
dt <- import("datetime")

#  Project settings --------
project_id <- "mpaeu"
taxons <- extract_storr()
taxons <- 100801
s3_path <- "s3://obis-maps/sdm"
sp_results_folder <- "results"
div_results_folder <- ""
hab_results_folder <- ""
catalog_output <- "stac"

# Root of the catalogue --------
root_catalog <- pystac$Catalog(
    id = "sdm-catalog",
    description = "OBIS species distribution modeling: maps for species, diversity and habitat."
)
root_catalog$describe()

# Themes -------
themes <- c("species", "diversity", "habitat")
theme_catalogs <- vector("list", length = length(themes))
names(theme_catalogs) <- themes

for (th in themes) {
    theme_cat <- pystac$Catalog(
        id = glue("{th}-catalog"),
        description = glue("{stringr::str_to_title(th)} modeling results")
    )
    root_catalog$add_child(theme_cat)
    theme_catalogs[[th]] <- theme_cat
}

# Project catalogue for species -----
species_project_catalog <- pystac$Catalog(
    id = glue("species-{project_id}"),
    description = glue("Species modeling results for model={project_id}")
)
theme_catalogs[["species"]]$add_child(species_project_catalog)

theme_catalogs$species$describe()

species_collection <- pystac$Collection(
    id = glue("species-{project_id}-collection"),
    description = glue("Species distribution models for model={project_id}"),
    extent = pystac$Extent(
        spatial = pystac$SpatialExtent(c(-180, -90, 180, 90)),
        temporal = pystac$TemporalExtent(list(
            dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),
            dt$datetime(2100L, 1L, 1L, tzinfo = dt$timezone$utc)
        ))
        # c(as.Date("2025-01-01"), as.Date("2100-01-01")))
    ),
    license = "CC-BY-4.0"
)
species_project_catalog$add_child(species_collection)
species_project_catalog$describe()

for (sp in taxons) {
    # PREDICTIONS ------
    predictions <- file.path(
        sp_results_folder,
        glue("taxonid={sp}/model={project_id}/predictions")
    ) |>
        list.files(pattern = "\\.tif$")
    predictions <- gsub("mask", "what=mask", predictions) #temp

    mask_file <- predictions[grepl("what=mask", predictions)]
    predictions <- predictions[!grepl("what=mask", predictions)]

    mess_file <- predictions[grepl("what=mess", predictions)]
    shape_file <- predictions[grepl("what=shape", predictions)]
    therm_file <- predictions[grepl("what=thermenvelope", predictions)]

    predictions <- predictions[!grepl("what=mess|what=shape|what=thermenvelope", predictions)]

    fit_pts <- file.path(
        sp_results_folder,
        glue("taxonid={sp}/model={project_id}/taxonid={sp}_model={project_id}_what=fitocc.parquet")
    ) |>
        arrow::read_parquet()
    fit_extent <- c(min(fit_pts$decimalLongitude), min(fit_pts$decimalLatitude),
                    max(fit_pts$decimalLongitude), max(fit_pts$decimalLatitude))
    fit_extent <- round(fit_extent, 2)

    log_file <- jsonlite::read_json(
        file.path(
            sp_results_folder,
            glue("taxonid={sp}/model={project_id}/taxonid={sp}_model={project_id}_what=log.json")
        ), simplifyVector = TRUE
    )

    if (length(as.vector(log_file$model_good)) > 1) {
        good_models <- c(as.vector(log_file$model_good), "ensemble")
    } else {
        good_models <- as.vector(log_file$model_good)
    }

    item <- pystac$Item(
        id = glue("taxonid={sp}"),
        geometry = NULL,
        bbox = c(-180, -90, 180, 90),
        datetime = dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),
        properties = list(
            taxonid = as.integer(sp),
            model = project_id,
            scientificName = as.vector(log_file$scientificName),
            group = as.vector(log_file$group),
            hab_depth = as.vector(log_file$hab_depth),
            fit_bbox = fit_extent,
            fit_n_points = as.vector(log_file$model_fit_points),
            methods = good_models
        )
    )

    for (pred in predictions) {
        method <- gsub("method=", "", sub(".*?(method=.*?)_scen.*", "\\1", pred))
        method_name <- dplyr::case_when(
            method == "maxent" ~ "Maxent",
            method == "xgboost" ~ "XGBoost",
            method == "rf_classification_ds" ~ "Random Forest - Down-sampled classification",
            method == "ensemble" ~ "Ensemble of models"
        )
        if (grepl("what=bootcv", pred)) {
            scenario <- gsub("scen=", "", sub(".*?(scen=.*?)_what.*", "\\1", pred))
        } else {
            scenario <- gsub("scen=", "", sub(".*?(scen=.*?)_cog.*", "\\1", pred))
        }

        pred_type <- gsub(glue("taxonid={sp}_model={project_id}_"), "", pred)
        if (grepl("bootcv", pred_type)) {
            pred_type <- gsub("bootcv_cog.tif", "uncertainty", pred_type)
        } else {
            pred_type <- gsub("_cog.tif", "_what=prediction", pred_type)
        }

        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/predictions/{pred}"))

        item$add_asset(
            pred_type,
            pystac$Asset(
                href = file_path,
                media_type = pystac$MediaType$COG,
                roles = list("data")
            )
        )
    }

    # METRICS ------
    metrics <- file.path(
        sp_results_folder,
        glue("taxonid={sp}/model={project_id}/metrics")
    ) |>
        list.files()

    for (met in metrics) {
        if (tools::file_ext(met) == "rds") {
            met_type <- "application/octet-stream"
            method <- NULL
            method_name <- NULL
        } else if (tools::file_ext(met) == "parquet") {
            met_type <- "application/vnd.apache.parquet" #pystac not working here
            if (grepl("method=", met)) {
                method <- gsub("method=", "", sub(".*?(method=.*?)_what.*", "\\1", met))
                method_name <- dplyr::case_when(
                    method == "maxent" ~ "Maxent",
                    method == "xgboost" ~ "XGBoost",
                    method == "rf_classification_ds" ~ "Random Forest - Down-sampled classification",
                    method == "ensemble" ~ "Ensemble of models"
                )
            } else {
                method <- NULL
                method_name <- NULL
            }
        }

        item_what_short <- gsub("what=", "", sub(".*?(what=.*?)\\..*", "\\1", met))
        item_explanation <- dplyr::case_when(
            item_what_short == "cvmetrics" ~ "Cross-validation metrics",
            item_what_short == "fullmetrics" ~ "Metrics on the full data (for reference only, use always the CV metrics)",
            item_what_short == "respcurves" ~ "Partial response curves",
            item_what_short == "varimportance" ~ "Variables importance",
            item_what_short == "biasmetrics" ~ "Data bias metrics (K-function and L-function)",
            item_what_short == "posteval_hyperniche" ~ "Post-evaluation: hyperniche",
            item_what_short == "posteval_niche" ~ "Post-evaluation: niche overlap",
            item_what_short == "thresholds" ~ "Model thresholds",
            .default = item_what_short
        )
        item_what <- gsub(
            paste0("\\.", tools::file_ext(met)), "",
            gsub(glue("taxonid={sp}_model={project_id}_"), "", met)
        )

        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/metrics/{met}"))
        item$add_asset(
            item_what,
            pystac$Asset(
                href = file_path,
                title = item_what_short,
                description = item_explanation,
                media_type = met_type,
                roles = list("data")
            )
        )
    }

    # MODELS ------
    models <- file.path(
        sp_results_folder,
        glue("taxonid={sp}/model={project_id}/models")
    ) |>
        list.files()

    for (mod in models) {
        mod_type <- "application/octet-stream"
        method <- gsub("method=", "", sub(".*?(method=.*?)_what.*", "\\1", mod))
        method_name <- dplyr::case_when(
            method == "maxent" ~ "Maxent",
            method == "xgboost" ~ "XGBoost",
            method == "rf_classification_ds" ~ "Random Forest - Down-sampled classification",
            method == "ensemble" ~ "Ensemble of models"
        )

        item_what <- gsub("\\.rds", "", gsub(glue("taxonid={sp}_model={project_id}_"), "", mod))
        item_explanation <- glue("RDS (R serialized format) file containing the model fit for {method}")

        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/metrics/{mod}"))
        item$add_asset(
            item_what,
            pystac$Asset(
                href = file_path,
                title = method_name,
                description = item_explanation,
                media_type = mod_type,
                roles = list("data")
            )
        )
    }

    # OTHER RESOURCES ------
    other_resources <- file.path(
        sp_results_folder,
        glue("taxonid={sp}/model={project_id}/")
    ) |>
        list.files()
    other_resources <- other_resources[!other_resources %in% c("metrics", "predictions", "models", "figures")]

    other_resources <- c(
        other_resources,
        mask_file,
        mess_file,
        shape_file,
        therm_file
    )

    for (other in other_resources) {
        asset_type <- gsub(paste0("\\.", tools::file_ext(other)), "", gsub(
            glue("taxonid={sp}_model={project_id}_"), "", other
        ))
        asset_type <- gsub("_cog", "", asset_type)

        asset_explanation <- dplyr::case_when(
            grepl("what=fitocc", other) ~ "Points used to fit the model",
            grepl("what=log", other) ~ "Log file containing all model details",
            grepl("what=mask", other) ~ "Mask to restrict predictions (multi band, each band is a different mask)",
            grepl("what=mess", other) ~ "MESS uncertainty metric",
            grepl("what=shape", other) ~ "SHAPE uncertainty metric",
            grepl("what=thermenvelope", other) ~ "Thermal envelope (each band is a time period)"
        )

        asset_file_type <- dplyr::case_when(
            grepl("what=fitocc", other) ~ "application/vnd.apache.parquet",
            grepl("what=log", other) ~ pystac$MediaType$JSON,
            grepl("what=mask", other) ~ pystac$MediaType$COG,
            grepl("what=mess", other) ~ pystac$MediaType$COG,
            grepl("what=shape", other) ~ pystac$MediaType$COG,
            grepl("what=thermenvelope", other) ~ pystac$MediaType$COG
        )

        file_sub_path <- dplyr::case_when( 
            grepl("what=fitocc", other) ~ glue("species/taxonid={sp}/model={project_id}/{other}"),
            grepl("what=log", other) ~ glue("species/taxonid={sp}/model={project_id}/{other}"),
            grepl("what=mask", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=mess", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=shape", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=thermenvelope", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}")
        )
        file_sub_path <- file.path(s3_path, file_sub_path)

        item$add_asset(
            asset_type,
            pystac$Asset(
                href = file_sub_path,
                media_type = asset_file_type,
                title = gsub("what=", "", asset_type),
                description = asset_explanation,
                roles = list("data")
            )
        )
    }

    species_collection$add_item(item)
}

#species_collection$describe()
#root_catalog$describe()

root_catalog$normalize_and_save(
    root_href = catalog_output,
    catalog_type = pystac$CatalogType$SELF_CONTAINED
)
