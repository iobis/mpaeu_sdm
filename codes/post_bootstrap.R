############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
##################### SDM - bootstrap results of models ########################

bootstrap_model <- function(species, iterations = 20,
                            model_acro = "mpaeu",
                            results_folder = "results",
                            target = "all",
                            n_back_max = NULL) {

    source("functions/components_model_species.R")
    log_file <- jsonlite::read_json(glue::glue("{results_folder}/taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=log.json"))

    best_conf <- log_file[["model_bestparams"]]
    best_variables <- unlist(log_file[["model_details"]][["variables"]])
    good_models <- unlist(log_file[["model_good"]])
    metric_type <- log_file[["model_good_metric"]][[1]]
    group <- log_file[["group"]][[1]]

    if (target == "all") {
        to_do <- good_models
    } else if (target == "best") {
        metrics <- lapply(good_models, function(mn) {
            mfiles <- list.files(glue::glue("results/taxonid={species}/model={model_acro}/metrics/"), full.names = T)
            mf <- mfiles[grepl("what=cvmetr", mfiles)]
            mf <- mf[grepl(paste0("method=", substr(mn, 1, 2)), mf)]
            mo <- arrow::read_parquet(mf)
            return(mean(mo[metric_type][, 1], na.rm = T))
        })
        to_do <- good_models[which.max(metrics)]
    } else {
        to_do <- base::intersect(target, good_models)
    }

    if (length(to_do) < 1) {
        cli::cli_abort("No good model corresponding to target='{target}'. Available: {.val {good_models}} model{?s}")
    }

    # For back compatibility
    if ("hab_depth" %in% names(log_file)) {
        hab_depth <- log_file[["hab_depth"]][[1]]
    } else {
        eco_info <- arrow::open_csv_dataset("data/species_ecoinfo.csv") %>%
            filter(taxonID == species) %>%
            collect()

        if (nrow(eco_info) < 1) {
            eco_info <- obissdm::mp_get_ecoinfo(
                species_list = species,
                outfile = NULL,
                return_table = TRUE,
                show_progress = FALSE
            )
        }

        # Select ecological information for the species
        hab <- eco_info$mode_life
        hab_depth <- hab_to_depth(hab)
    }

    # Load env layers
    env <- obissdm::get_envofgroup(group,
        depth = hab_depth, load_all = T,
        conf_file = "sdm_conf.yml", verbose = FALSE
    )
    env <- terra::subset(env$layers, best_variables)

    mask_l <- terra::rast(glue::glue("results/taxonid={species}/model={model_acro}/predictions/taxonid={species}_model={model_acro}_mask_cog.tif"))
    mask_l <- terra::subset(mask_l, "fit_region")
    mask_l[mask_l == 0] <- NA

    env_masked <- terra::mask(env, mask_l)

    # Load points
    pts <- arrow::read_parquet(glue::glue("results/taxonid={species}/model={model_acro}/taxonid={species}_model={model_acro}_what=fitocc.parquet"))

    quad_n <- log_file[["model_details"]][["background_size"]][[1]]

    if (!is.null(n_back_max)) {
        quad_n <- ifelse(quad_n > n_back_max, n_back_max, quad_n)
    }

    env_xy <- terra::as.data.frame(prod(env_masked), xy = T)[,1:2]
    env_xy_samp <- sample(1:nrow(env_xy), size = quad_n)
    env_xy <- env_xy[env_xy_samp,]
    colnames(env_xy) <- c("decimalLongitude", "decimalLatitude")

    sampled_quad <- cbind(
        env_xy, presence = 0,
        terra::extract(env_masked, env_xy, ID = F)
    )

    for (i in 1:iterations) {

        samp_pts <- sample(1:nrow(pts), nrow(pts), replace = T)
        samp_quad <- sample(1:nrow(sampled_quad), nrow(sampled_quad), replace = T)

        sp_data <- obissdm::mp_prepare_data(pts[samp_pts,], eval_data = NULL,
                               species_id = unlist(log_file[["scientificName"]]),
                               env_layers = env_masked,
                               pred_quad = sampled_quad[samp_quad,],
                               verbose = FALSE)

        for (z in 1:length(to_do)) {
            fit <- obissdm::model_bootstrap()
        }

    }





}

listviewer::jsonedit(log_file)

# END
