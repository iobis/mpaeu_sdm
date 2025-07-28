limit_by_closest <- function(
    results_folder,
    model_acro,
    species,
    algo,
    scenario,
    max_dispersal = 100000,
    target_th = "p10",
    target_mask = "fit_region_max_depth",
    verbose = TRUE
) {

    require(dplyr)
    if (verbose) cli::cli_progress_step("Processing species {.val {species}}, model {.val {algo}}, scenario {.val {scenario}}, with maximum dispersal {.val {max_dispersal/1000}} km")

    # Load prediction
    pred <- terra::rast(file.path(
        results_folder, glue::glue("taxonid={species}/model={model_acro}/predictions"),
        glue::glue("taxonid={species}_model={model_acro}_method={algo}_scen={scenario}_cog.tif")
    ))

    # Load points and apply buffer
    fit_pts <- arrow::read_parquet(file.path(
        results_folder,
        glue::glue("taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=fitocc.parquet")
    ))
    fit_pts <- terra::vect(as.data.frame(fit_pts), geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    fit_pts_buffer <- terra::buffer(fit_pts, max_dispersal)

    # Load masks
    masks <- terra::rast(file.path(
        results_folder, glue::glue(
            "taxonid={species}/model={model_acro}/predictions/taxonid={species}_model={model_acro}_mask_cog.tif"
        )
    ))
    masks <- subset(masks, target_mask)
    terra::NAflag(masks) <- 0

    # Apply mask
    pred <- terra::mask(pred, masks)

    # Apply threshold
    thresholds <- arrow::read_parquet(file.path(
        results_folder,
        glue::glue(
            "taxonid={species}/model={model_acro}/metrics/taxonid={species}_model={model_acro}_what=thresholds.parquet"
        )
    ))
    th <- thresholds |>
        filter(model == ifelse(grepl("rf_", algo), "rf", algo)) |>
        pull(target_th) |>
        (\(x) {x * 100})()
    pred_th <- terra::classify(pred, matrix(c(-Inf, th, NA), ncol = 3), right = F)

    # Mask 0s
    pred_masked <- terra::mask(pred_th, pred_th != 0, maskvalues = 0)

    # Label patches and get polygon
    pred_patches <- terra::patches(pred_masked, directions = 8)
    patches_pols <- terra::as.polygons(pred_patches, dissolve = TRUE)

    # Find intersections
    touched_pt <- terra::intersect(patches_pols, fit_pts_buffer)
    touched_pt_ids <- unique(terra::values(touched_pt)[[1]])

    # Mask original
    touched_pt_f <- terra::classify(pred_patches, rcl = cbind(setdiff(unique(terra::values(pred_patches)), touched_pt_ids), NA))
    touched_pt_f <- terra::classify(touched_pt_f, matrix(c(-Inf, Inf, 1), ncol = 3))
    touched_pt_f <- terra::as.int(touched_pt_f)

    if (verbose) cli::cli_progress_done()
    return(touched_pt_f)

}

check_good_models <- function(sp, model_acro, results_folder, target_models) {
    good_models <- jsonlite::read_json(
        file.path(results_folder,
                  glue::glue("taxonid={sp}/model={model_acro}/taxonid={sp}_model={model_acro}_what=log.json"))
    )
    good_models <- unlist(good_models$model_good, use.names = F)
    if (is.numeric(good_models)) {
        good_models <- "esm"
    }
    good_models <- gsub("rf", "rf_ds_classification", good_models)
    intersect(good_models, target_models)
}
