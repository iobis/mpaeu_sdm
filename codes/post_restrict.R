############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################## Post-processing #################################
#################### Create dispersal restricting masks ########################
source("functions/post_proc_functions.R")
verbosity <- TRUE

results_folder <- "/data/scps/v5/results"
model_acro <- "mpaeu"
species <- 124287
models <- c("maxent", "rf_ds_classification", "xgboost")
scenarios <- c("current", 
               paste0(rep(c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"), 2),
                      rep(c("_dec50", "_dec100"), each = 5)))

for (sp in species) {
    good_models <- check_good_models(sp, model_acro, results_folder, models)
    if (length(good_models) < 1) next
    for (m in good_models) {
        out_rast <- vector(mode = "list", length = length(scenarios))
        for (s in seq_along(scenarios)) {
            out_rast[[s]] <- try(limit_by_closest(
                results_folder, model_acro, species = sp,
                algo = m, scenario = scenarios[s], target_mask = "fit_region",
                verbose = verbosity
            ))
        }
        out_rast <- rast(out_rast)
        names(out_rast) <- scenarios
        outf <- file.path(results_folder,
                          glue::glue("taxonid={sp}/model={model_acro}/predictions"),
                          glue::glue("taxonid={sp}_model={model_acro}_method={m}_what=postrestrict.tif"))
        writeRaster(as.int(out_rast), outf, datatype = "INT1U", overwrite = T)
        obissdm::cogeo_optim(outf, verbose = verbosity)
    }
}
# END
