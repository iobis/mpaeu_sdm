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
s3_path <- "s3//obis-maps/sdm"
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

        pred_type <- ifelse(grepl("what=bootcv", pred), "uncertainty", "prediction")
        
        item_id <- glue("taxonid={sp}_model={project_id}_method={method}_scen={scenario}_what={pred_type}")
        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/predictions/{pred}"))

        item <- pystac$Item(
            id = item_id,
            geometry = NULL,
            bbox = c(-180, -90, 180, 90),
            datetime = dplyr::case_when(
                scenario == "current" ~ dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),#as.Date("2025-01-01"),
                grepl("dec50", scenario) ~ dt$datetime(2050L, 1L, 1L, tzinfo = dt$timezone$utc),#as.Date("2050-01-01"),
                grepl("dec100", scenario) ~ dt$datetime(2100L, 1L, 1L, tzinfo = dt$timezone$utc)#as.Date("2100-01-01")
            ),
            properties = list(
                taxonid = sp,
                model = project_id,
                method = method,
                method_full_name = method_name,
                scenario = scenario,
                fit_bbox = fit_extent
            )
        )

        item$add_asset(
            pred_type,
            pystac$Asset(
                href = file_path,
                media_type = pystac$MediaType$COG,
                roles = list("data")
            )
        )

        species_collection$add_item(item)
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
        

        item_what <- gsub("what=", "", sub(".*?(what=.*?)\\..*", "\\1", met))
        item_explanation <- dplyr::case_when(
            item_what == "cvmetrics" ~ "Cross-validation metrics",
            item_what == "fullmetrics" ~ "Metrics on the full data (for reference only, use always the CV metrics)",
            item_what == "respcurves" ~ "Partial response curves",
            item_what == "varimportance" ~ "Variables importance",
            item_what == "biasmetrics" ~ "Data bias metrics (K-function and L-function)",
            item_what == "posteval_hyperniche" ~ "Post-evaluation: hyperniche",
            item_what == "posteval_niche" ~ "Post-evaluation: niche overlap",
            item_what == "thresholds" ~ "Model thresholds",
            .default = item_what
        )

        item_id <- gsub(paste0("\\.", tools::file_ext(met)), "", basename(met))
        item <- pystac$Item(
            id = item_id,
            geometry = NULL,
            bbox = NULL,
            datetime = dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),#as.Date("2025-01-01"),
            properties = list(
                taxonid = sp,
                model = project_id,
                method = method,
                method_full_name = method_name,
                metric_type = item_what,
                metric_explanation = item_explanation
            )
        )

        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/metrics/{met}"))
        item$add_asset(
            "metric",
            pystac$Asset(
                href = file_path,
                media_type = met_type,
                roles = list("data")
            )
        )

        species_collection$add_item(item)
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

        item_what <- "model"
        item_explanation <- glue("RDS (R serialized format) file containing the model fit for {method}")

        item_id <- gsub(paste0("\\.", tools::file_ext(mod)), "", basename(mod))
        item <- pystac$Item(
            id = item_id,
            geometry = NULL,
            bbox = NULL,
            datetime = dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),#as.Date("2025-01-01"),
            properties = list(
                taxonid = sp,
                model = project_id,
                method = method,
                method_full_name = method_name,
                metric_type = item_what,
                metric_explanation = item_explanation
            )
        )

        file_path <- file.path(s3_path, glue("species/taxonid={sp}/model={project_id}/metrics/{mod}"))
        item$add_asset(
            "metric",
            pystac$Asset(
                href = file_path,
                media_type = mod_type,
                roles = list("data")
            )
        )

        species_collection$add_item(item)
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
        asset_type <- dplyr::case_when(
            grepl("what=fitocc", other) ~ "fit_points",
            grepl("what=log", other) ~ "log_file",
            grepl("what=mask", other) ~ "predictions_mask",
            grepl("what=mess", other) ~ "mess_uncertainty",
            grepl("what=shape", other) ~ "shape_uncertainty",
            grepl("what=thermenvelope", other) ~ "thermal_envelope"
        )
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
        asset_bbox <- dplyr::case_when(
            grepl("what=fitocc", other) ~ fit_extent,
            grepl("what=log", other) ~ NA,
            .default = c(-180, -90, 180, 90)
        )
        if (is.na(asset_bbox[1])) asset_bbox <- NULL

        file_sub_path <- dplyr::case_when( 
            grepl("what=fitocc", other) ~ glue("species/taxonid={sp}/model={project_id}/{other}"),
            grepl("what=log", other) ~ glue("species/taxonid={sp}/model={project_id}/{other}"),
            grepl("what=mask", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=mess", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=shape", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}"),
            grepl("what=thermenvelope", other) ~ glue("species/taxonid={sp}/model={project_id}/predictions/{other}")
        )
        file_sub_path <- file.path(s3_path, file_sub_path)

        item_id <- gsub(
            "_cog", "", gsub(paste0("\\.", tools::file_ext(other)), "", basename(other))
        )
        item <- pystac$Item(
            id = item_id,
            geometry = NULL,
            bbox = asset_bbox,
            datetime = dt$datetime(2025L, 1L, 1L, tzinfo = dt$timezone$utc),#as.Date("2025-01-01"),
            properties = list(
                taxonid = sp,
                model = project_id,
                object_type = asset_type,
                object_explanation = asset_explanation
            )
        )

        item$add_asset(
            asset_type,
            pystac$Asset(
                href = file_sub_path,
                media_type = asset_file_type,
                roles = list("data")
            )
        )

        species_collection$add_item(item)
    }
}

#species_collection$describe()
#root_catalog$describe()

root_catalog$normalize_and_save(
    root_href = catalog_output,
    catalog_type = pystac$CatalogType$SELF_CONTAINED
)



# Project (MPA Europe)
project_catalog <- pystac$Catalog(
  id = "project-mpaeu",
  description = "Species distribution models for the MPA Europe project"
)
root_catalog$add_child(project_catalog)
root_catalog$describe()

Extent <- pystac$Extent
SpatialExtent <- pystac$SpatialExtent
TemporalExtent <- pystac$TemporalExtent

species_collection <- pystac$Collection(
  id = "taxonid=12345",
  description = "Species distribution models for taxon 12345",
  extent = Extent(
    spatial = SpatialExtent(list(list(-180, -90, 180, 90))),
    temporal = TemporalExtent(list(list(NULL, NULL)))
  ),
  license = "CC-BY-4.0"
)
project_catalog$add_child(species_collection)

Asset <- pystac$Asset
MediaType <- pystac$MediaType

species_collection$add_asset(
  "mask",
  Asset(
    href = "s3://yourbucket/results/taxonid=12345/mask.tif",
    media_type = MediaType$COG,
    roles = list("mask")
  )
)
Item <- pystac$Item

item <- Item(
  id = "taxonid=12345_rf_current",
  geometry = NULL,
  bbox = list(-180, -90, 180, 90),
  datetime = as.POSIXct(Sys.Date()),
  properties = dict(
    taxonid = "12345",
    model = "mpaeu",
    method = "rf",
    scenario = "current"
  )
)

item$add_asset(
  "prediction",
  Asset(
    href = "s3://yourbucket/results/taxonid=12345/model=mpaeu/taxonid=12345_model=mpaeu_method=rf_scen=current_cog.tif",
    media_type = MediaType$COG,
    roles = list("data")
  )
)

species_collection$add_item(item)

CatalogType <- pystac$CatalogType
root_catalog$normalize_and_save(
  root_href = "stac_r_reticulate",
  catalog_type = CatalogType$SELF_CONTAINED
)


library(rstac)
s_obj <- stac("stac_r_reticulate/catalog.json", force_version = T)
get_request(s_obj)


root_catalog = pystac$Catalog$from_file("stac_r_reticulate/catalog.json")
root_catalog$describe()
