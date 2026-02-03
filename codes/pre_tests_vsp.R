############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ SDM - testing methods #############################

# Load packages ----
library(obissdm)
library(arrow)
library(terra)
library(future)
library(furrr)
library(dplyr)
library(storr)
source("functions/load_env.R")
source("functions/components_model_species.R")
source("functions/auxiliary_modelfit.R")
n_workers <- 4
plan(multisession, workers = n_workers)
set.seed(2023)

# Prepare settings ----
st <- storr_rds("vsptest_storr")
outfolder <- "results/vsp_testing/"
fs::dir_create(outfolder)

# Prepare functions -------
# Prepare fit function for parallel processing
fit_parallel <- function(species, variables, storr_control, n_workers, verbose = TRUE) {
    terra::terraOptions(memfrac = (0.7 / n_workers))

    files <- list.files(file.path("data/virtual_species", paste0("key=", species)),
        full.names = T
    )
    files <- files[grepl("occurrences", files)]

    modes <- unique(gsub("_rep*.*", "", basename(files)))

    if (verbose) cli::cli_alert_info("{length(modes)} data types with {sum(grepl(modes[1], files))/2} replicates.")

    for (mo in modes) {
        if (verbose) cli::cli_alert_info("Executing mode {mo}")

        files_sel <- files[grepl(mo, files)]
        files_sel_pa <- files_sel[grepl("_pa", files_sel)]
        files_sel <- files_sel[!grepl("_pa", files_sel)]

        for (tf in seq_along(files_sel)) {
            tf_base <- files_sel[tf]
            tf_pa <- files_sel_pa[tf]

            run_code <- paste0(species, "_", mo, "_", tf)

            status <- try(fit_model_vsp(run_code, tf_base, tf_pa, variables, outfolder, verbose = verbose),
             silent = FALSE)

            if (inherits(status, "try-error")) {
                storr_control$set(run_code, "failed")
            } else {
                storr_control$set(run_code, "success")
            }
        }
    }
    return(invisible(NULL))
}

# Prepare model fit function, with same parameters used in the model fitting with real species
fit_model_vsp <- function(species, base_file, pa_file, variables, outfolder, verbose) {

    outacro <- "vsp"

    # Algorithms to be used
    algos <- c("maxent", "rf", "xgboost")
    # Personalized options
    algo_opts <- obissdm::sdm_options()[algos]
    algo_opts$maxent$features <- c("lq", "h")
    algo_opts$maxent$remult <- seq_len(4)
    algo_opts$xgboost$gamma <- c(0, 4)
    algo_opts$xgboost$shrinkage <- c(0.1, 0.3)
    algo_opts$xgboost$scale_pos_weight <- c("balanced", "equal")

    # Should areas be masked by the species depth?
    limit_by_depth <- TRUE
    # A buffer to be applied to the depth limitation
    depth_buffer <- 50
    # Assess spatial bias?
    assess_bias <- FALSE
    # Quadrature size
    quad_samp <- 0.01 # 1% of the total number of points
    # Target metric
    tg_metric <- "cbi"
    # Metric threshold
    tg_threshold <- 0.3

    post_eval <- c("sst")

    if (verbose) {
        verb_1 <- verb_2 <- TRUE
    } else {
        verb_1 <- verb_2 <- FALSE
    }

    if (verb_1) cli::cli_alert_info("Starting model for species {species}")

    # Record timing
    treg <- obissdm::.get_time()

    # Define global objects and fill log
    coord_names <- c("decimalLongitude", "decimalLatitude")
    model_log <- obissdm::gen_log(algos)
    model_log$taxonID <- species
    model_log$group <- NULL
    model_log$model_acro <- "vsp"
    model_log$model_date <- Sys.Date()
    if (!is.null(algo_opts)) {
        algo_opts_nams <- names(algo_opts)
        for (opt in algo_opts_nams) model_log$algorithms_parameters[[opt]] <- algo_opts[[opt]]
    }
    algorithms <- algos

    # Define output paths
    metric_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
    pred_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/")
    mod_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
    fig_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/figures/")

    if (verb_1) cli::cli_alert_info("Reading data")
    species_data <- arrow::read_parquet(base_file)
    fit_pts <- species_data[, 1:2] |> as.data.frame()

    eval_pts <- arrow::read_parquet(pa_file)[, 1:3] |> as.data.frame()

    treg <- obissdm::.get_time(treg, "Species data loading")

    model_log$scientificName <- species

    model_log$hab_depth <- "surface"
    hab_depth <- "surf"

    # Load environmental data
    env <- load_env(variables, terrain_vars = c("bathymetry_mean"))

    env <- list(hypothesis = list(base = names(env)), layers = env)

    env <- .cm_check_coastal(species_data, env, coord_names, verb_2)


    # Load Realms shapefile (we call ecoregions for convenience)
    ecoregions <- vect("data/shapefiles/MarineRealms_BO.shp")

    model_log$n_init_points <- nrow(fit_pts)

    treg <- obissdm::.get_time(treg, "Data loading")



    # PART 2: DATA PREPARING ----
    if (verb_1) cli::cli_alert_info("Preparing data")

    # Check which eco-regions are covered by the points
    prep_eco <- .cm_check_ecoregions(
        ecoregions, fit_pts, eval_pts[,1:2], env, limit_by_depth, depth_buffer,
        coord_names, model_log, verb_1
    )
    env <- prep_eco$env
    model_log <- prep_eco$model_log
    bath_pts <- prep_eco$bath_pts
    ecoreg_occ <- prep_eco$ecoreg_occ
    ecoreg_sel <- prep_eco$ecoreg_sel
    max_depth <- min(bath_pts)
    rm(prep_eco)

    # Check if species is from shallow/coastal areas and remove distcoast/bathymetry
    env <- .cm_check_shallow(bath_pts, env, verb_1)

    # Assess SAC
    fit_pts_sac <- try(obissdm::outqc_sac_mantel(fit_pts,
        env_layers = terra::subset(env$layers, env$hypothesis[[1]]),
        plot_result = FALSE,
        verbose = verb_2
    ))

    if (inherits(fit_pts_sac, "try-error")) {
        fit_pts_sac <- fit_pts
    }

    # Make data object
    quad_n <- .cm_calc_quad(env, quad_samp, fit_pts_sac)
    model_log$model_details$background_size <- quad_n

    sp_data <- mp_prepare_data(fit_pts_sac,
        eval_data = eval_pts,
        species_id = species,
        env_layers = env$layers,
        quad_number = quad_n,
        verbose = verb_2
    )

    block_grid <- get_block_grid(sp_data, env$layers,
        sel_vars = env$hypothesis$basevars,
        verbose = verb_2
    )

    sp_data <- mp_prepare_blocks(sp_data,
        method = "manual",
        block_types = "spatial_grid",
        n_iterate = 300,
        retry_if_zero = TRUE,
        manual_shp = block_grid,
        verbose = verb_2
    )

    if (any(table(sp_data$training$presence, sp_data$blocks$folds[["spatial_grid"]])[2, ] == 0)) {
        stop("Blocks with less than 1 point. Failed.")
    }

    model_log$model_details$block_size <- unname(sp_data$blocks$grid_resolution)
    model_log$model_fit_points <- sum(sp_data$training$presence)
    model_log$model_eval_points <- sum(sp_data$eval_data$presence)

    treg <- obissdm::.get_time(treg, "Data preparing")

    # PART 3: ASSESS SPATIAL BIAS ----
    if (assess_bias) {
        if (verb_1) cli::cli_alert_info("Assessing spatial bias")

        .cm_spat_bias(metric_out, species, outacro, env, sp_data)

        treg <- obissdm::.get_time(treg, "Sample bias assessment")
    }

    # PART 4: MODEL SELECTION ----
    if (verb_1) cli::cli_alert_info("Model selection")

    if (length(env$hypothesis) != 1) {
        multi_mod <- sdm_multhypo(
            sdm_data = sp_data, sdm_method = "rf",
            variables = env$hypothesis,
            options = algo_opts[["rf"]],
            verbose = verb_2
        )
    } else {
        if (verb_1) cli::cli_alert_info("Only one hypothesis available - skipping")
        multi_mod <- list(
            best_model = names(env$hypothesis),
            best_variables = env$hypothesis[[1]]
        )
    }

    model_log$model_details$hypothesis_tested <- env$hypothesis
    model_log$model_details$best_hypothesis <- multi_mod$best_model
    model_log$model_details$variables <- multi_mod$best_variables

    # Prepare data for the best model
    sp_data$training <- sp_data$training[, c("presence", multi_mod$best_variables)]
    sp_data$eval_data <- sp_data$eval_data[, c("presence", multi_mod$best_variables)]

    treg <- obissdm::.get_time(treg, "Model selection")


    # PART 5: FIT MODELS ----
    if (verb_1) cli::cli_alert_info("Starting model fitting")
    fit_obj <- .cm_model_fit(
        algorithms, algo_opts, sp_data,
        model_log, verb_1, verb_2
    )
    model_fits <- fit_obj$fits
    model_log <- fit_obj$logs
    rm(fit_obj)
    treg <- obissdm::.get_time(treg, "Model fit")

    if (any(unlist(model_log$model_result) == "succeeded")) {
        if (verb_1) cli::cli_alert_info("Model fitting concluded. Checking if models are good")

        # PART 6: PREDICTION ----
        good_models <- .cm_check_good_models(model_fits, tg_metric, tg_threshold)

        if (length(good_models) > 0) {
            if (verb_1) cli::cli_alert_info("Good models available, predicting.")
            # Predict models
            model_predictions <- .cm_predict_models_vsp(
                good_models, model_fits,
                multi_mod, pred_out, group,
                hab_depth, sp_data, outfolder,
                outacro, species, env, variables,
                verb_1, verb_2
            )

            treg <- obissdm::.get_time(treg, "Model prediction")

            model_log$model_good <- algorithms[good_models]
            model_log$model_good_metric <- tg_metric
            model_log$model_good_threshold <- tg_threshold



            # PART 7: GET RESPONSE CURVES ----
            if (verb_1) cli::cli_alert_info("Getting response curves")
            # Predict models
            model_respcurves <- .cm_get_respcurves(
                good_models, model_fits, sp_data,
                env, multi_mod, metric_out, fig_out,
                species, outacro
            )

            treg <- obissdm::.get_time(treg, "Response curves")


            # Variables importance
            model_varimport <- .cm_get_importance(
                good_models, model_fits,
                sp_data, metric_out,
                species, outacro
            )

            treg <- obissdm::.get_time(treg, "Variable importance")

            # PART 8: ENSEMBLE ----
            if (length(good_models) > 1) {
                ens_obj <- .cm_ensemble_models(
                    species, good_models, algorithms, outfolder,
                    sp_data, model_fits, metric_out,
                    pred_out, outacro, model_log, model_predictions, verb_1
                )
                model_predictions <- ens_obj$preds
                model_log <- ens_obj$logs
                rm(ens_obj)
                treg <- obissdm::.get_time(treg, "Ensemble")
            }

            # PART 9: SAVE BINARIZATION INFO ----
            if (verb_1) cli::cli_alert_info("Saving binarization info")
            thresh_p10_mtp <- .cm_save_bin_info(
                model_predictions, sp_data,
                algorithms, good_models,
                metric_out, species, outacro
            )


            # PART 10: CREATE MASKS AND SAVE ----
            if (verb_1) cli::cli_alert_info("Saving masks")
            .cm_save_masks(
                ecoreg_occ, ecoreg_sel, multi_mod,
                env, model_predictions, sp_data,
                pred_out, species, outacro, coord_names, max_depth
            )
            treg <- obissdm::.get_time(treg, "Masks")



            # PART 11: EVALUATE MODELS USING OTHER TECHNIQUES ----
            if (verb_1) cli::cli_alert_info("Performing post-evaluation")
            post_ob <- .cm_check_posteval_n(sp_data, model_log)
            model_log <- post_ob$logs
            sp_data_post <- post_ob$dat
            rm(post_ob)


            model_log <- .cm_posteval_all(
                post_eval, model_predictions, sp_data_post,
                thresh_p10_mtp, algorithms, hab_depth,
                good_models, model_log, metric_out,
                species, outacro,
                model_varimport, env, verb_1
            )
            treg <- obissdm::.get_time(treg, "Post-evaluation")


            # PART 12: SAVE MODELS OBJECTS ----
            if (verb_1) cli::cli_alert_info("Saving models objects")

            # Save metrics
            .cm_save_metrics(
                model_fits, metric_out, species,
                outacro, good_models
            )

            # Update log with best parameters
            for (al in seq_along(algorithms)) {
                if (!is.logical(model_fits[[al]])) {
                    model_log$model_bestparams[[al]] <- model_fits[[al]]$parameters
                }
            }

            # Save models
            .cm_save_models(model_fits, outfolder, species, outacro, good_models)

            # Save fit points
            arrow::write_parquet(
                sp_data$coord_training[sp_data$training$presence == 1, ],
                paste0(
                    outfolder, "/taxonid=", species, "/model=", outacro,
                    "/taxonid=", species, "_model=", outacro, "_",
                    "what=fitocc.parquet"
                )
            )

            treg <- as.data.frame(unclass(treg))
            treg <- data.frame(what = row.names(treg), time_mins = treg[, 1])
            model_log$timings <- treg

            outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/")
            outfile <- paste0(
                outpath,
                "taxonid=", species, "_model=", outacro, "_what=log.json"
            )
            save_log(model_log, outfile)

            cli::cli_alert_success("Model concluded for species {species}.")

            st_status <- "succeeded"
        } else {
            if (verb_1) cli::cli_alert_warning("No good model available")
            st_status <- "no_good_model"
            if (dir.exists(file.path(outfolder, paste0("taxonid=", species)))) {
                fs::dir_delete(file.path(outfolder, paste0("taxonid=", species)))
            }
        }
    } else {
        if (verb_1) cli::cli_alert_warning("All models failed")
        st_status <- "all_failed"
        if (dir.exists(file.path(outfolder, paste0("taxonid=", species)))) {
            fs::dir_delete(file.path(outfolder, paste0("taxonid=", species)))
        }
    }
    return(invisible(NULL))
}

# Adapted "prediction function" adjusted for efficiency for the purpose of VSP testing
.cm_predict_models_vsp <- function(good_models, model_fits, multi_mod,
                                   pred_out, group, hab_depth,
                                   sp_data, outfolder, outacro, species, env, variables,
                                   verb_1, verb_2) {
    lapply(good_models, function(id) {
        model_name <- model_fits[[id]]$name
        if (verb_1) cli::cli_alert_info("Predicting model {id} - {model_name}")

        best_hyp <- multi_mod$best_model

        if (!dir.exists(pred_out)) fs::dir_create(pred_out)

        outmess <- paste0(pred_out, "taxonid=", species, "_model=", outacro, "_what=mess.tif")
        outshape <- paste0(pred_out, "taxonid=", species, "_model=", outacro, "_what=shape.tif")

        if (!file.exists(gsub("mess", "mess_cog", outmess))) {
            do_shape <- do_mess <- TRUE
        } else {
            do_shape <- do_mess <- FALSE
        }

        scenarios <- data.frame( # Here we predict to only part of the scenarios, what speed up the process
            scenario = c("current", rep(c("ssp126", "ssp585"),
                each = 1
            )),
            year = c(NA, rep(c("dec100"), 2))
        )

        for (k in seq_len(nrow(scenarios))) {
            if (is.na(scenarios$year[k])) {
                period <- NULL
            } else {
                period <- scenarios$year[k]
            }

            if (verb_1) cli::cli_alert_info("Predicting scenario {k} of {nrow(scenarios)}.")
            outpred <- paste0(
                pred_out, "taxonid=", species, "_model=", outacro,
                "_method=", model_name, "_scen=", scenarios$scenario[k],
                ifelse(is.null(period), "", paste0("_", period)), ".tif"
            )

            if (scenarios$scenario[k] == "current") {
                env_to_pred <- load_env(variables, terrain_vars = c("bathymetry_mean"))
                env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
            } else if (scenarios$scenario[k] == "ssp126") {
                env_to_pred <- load_env(variables, scenario = "ssp126", terrain_vars = c("bathymetry_mean"))
                env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
            } else {
                env_to_pred <- load_env(variables, scenario = "ssp585", terrain_vars = c("bathymetry_mean"))
                env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
            }

            pred <- predict(model_fits[[id]], env_to_pred)

            names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))

            pred <- pred * 100
            pred <- terra::as.int(pred)
            terra::writeRaster(pred, outpred, overwrite = T, datatype = "INT1U")
            cogeo_optim(outpred)

            if (k == 1) {
                pred_f <- pred
            }

            # Deactivate MESS/SHAPE to speed up
            do_mess <- do_shape <- FALSE

            if (do_mess) {
                # Save MESS
                if (verb_1) cli::cli_alert_info("Generating MESS map.")
                if (best_hyp == "coastal") {
                    to_mess <- terra::aggregate(env_to_pred, 12, na.rm = T)
                } else {
                    to_mess <- terra::aggregate(env_to_pred, 12)
                }
                mess_map <- ecospat::ecospat.mess(
                    na.omit(as.data.frame(to_mess, xy = T)),
                    cbind(sp_data$coord_training, sp_data$training[, 2:ncol(sp_data$training)])
                )
                mess_map_t <- to_mess[[1]]
                mess_map_t[] <- NA
                mess_map_t[terra::cellFromXY(mess_map_t, mess_map[, 1:2])] <- mess_map[, 5]
                mess_map <- mess_map_t

                names(mess_map) <- names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))

                if (k == 1) {
                    pred_mess <- mess_map
                } else {
                    pred_mess <- c(pred_mess, mess_map)
                }
            }

            if (do_shape) {
                if (verb_1) cli::cli_alert_info("Generating SHAPE map.")
                # Reduce dataset for faster implementing
                shape_data <- sp_data$training
                which_p <- which(shape_data$presence == 1)
                if (sum(shape_data$presence) > 1000) {
                    which_p <- sample(which_p, 1000)
                }
                which_a <- sample(which(shape_data$presence == 0), 1000)
                shape_data <- shape_data[c(which_p, which_a), ]
                shape_data_coords <- sp_data$coord_training[c(which_p, which_a), ]
                names(shape_data_coords) <- c("x", "y")

                shape_res <- try(flexsdm::extra_eval(
                    shape_data,
                    "presence",
                    projection_data = env_to_pred,
                    aggreg_factor = 12 # For faster implementing
                ), silent = T)

                if (!inherits(shape_res, "try-error")) {
                    if (k == 1) {
                        outpath_fig <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/figures/")
                        if (!dir.exists(outpath_fig)) fs::dir_create(outpath_fig)
                        outfile_fig <- paste0(
                            outpath_fig,
                            "taxonid=", species, "_model=", outacro, "_method=", model_name,
                            "_what=shape.png"
                        )
                        shape_plot <- suppressMessages(
                            try(flexsdm::p_extra(
                                training_data = cbind(shape_data_coords, shape_data),
                                pr_ab = "presence",
                                extra_suit_data = terra::aggregate(shape_res, 12),
                                projection_data = terra::aggregate(env_to_pred, 12),
                                geo_space = TRUE,
                                prop_points = 0.05,
                                alpha_p = 0.2
                            ), silent = T)
                        )

                        ragg::agg_png(outfile_fig, width = 6, height = 2.5, units = "in", res = 300, scaling = 0.4)
                        print(shape_plot)
                        dof <- dev.off()
                        rm(dof, shape_plot)
                    }

                    shape_res <- aggregate(shape_res, fact = 12)
                    shape_res <- terra::as.int(shape_res)
                    names(shape_res) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))

                    if (k == 1) {
                        pred_shape <- shape_res
                    } else {
                        pred_shape <- c(pred_shape, shape_res)
                    }
                }
            }
        }

        if (do_mess) {
            pred_mess <- terra::as.int(pred_mess)
            terra::writeRaster(pred_mess, outmess, overwrite = T, datatype = "INT1U")
            cogeo_optim(outmess)
        }
        if (do_shape) {
            terra::writeRaster(pred_shape, outshape, overwrite = T, datatype = "INT2U")
            cogeo_optim(outshape)
        }

        return(pred_f)
    })
}

# Apply models -----
env_vars <- c("thetao-mean", "sws-mean", "so-mean", "o2-mean")

result <- future_map(c(1001, 1002, 1003, 1004),
    fit_parallel,
    variables = env_vars, verbose = FALSE,
    storr_control = st, n_workers = n_workers,
    .progress = TRUE, .options = furrr_options(seed = T)
)

fit_parallel(
    1001,
    variables = env_vars, verbose = TRUE,
    storr_control = st, n_workers = n_workers)

# Check overlap -----
check_overlap <- function(species, output_folder = "results/vsp_testing") {
    results <- list()
    for (i in seq_along(species)) {
        spf <- list.files(output_folder)
        spf <- spf[grepl(species[i], spf)]
        
        sub_results <- list()
        for (k in seq_along(spf)) {
            models <- c("maxent", "rf", "xgboost")

            models_results <- list()
            for (mm in seq_along(models)) {
                sel_model <- models[mm]
                if (sel_model == "rf") {
                    sel_model <- "rf_classification_ds"
                }
                r1 <- rast(paste0(
                    "data/virtual_species/key=", species[i], "/suitability_key", species[i], ".tif"
                ))
                if (!file.exists(paste0(
                    "results/vsp_testing/taxonid=", species[i], "_occurrences_high_bias_1/model=vsp/predictions/taxonid=", species[i], "_occurrences_high_bias_1_model=vsp_method=", sel_model, "_scen=current_cog.tif"
                ))) next
                r2 <- rast(paste0(
                    "results/vsp_testing/taxonid=", species[i], "_occurrences_high_bias_1/model=vsp/predictions/taxonid=", species[i], "_occurrences_high_bias_1_model=vsp_method=", sel_model, "_scen=current_cog.tif"
                ))

                # thresh <- arrow::read_parquet(
                #     paste0("results/vsp_testing/", spf[k], "/model=vsp/metrics/", spf[k], "_model=vsp_what=thresholds.parquet")
                # )
                # thresh <- thresh |> filter(model == sel_model) |> select(p10) |> pull()

                #r2[r2 < (thresh * 100)] <- 0
                if (!file.exists(paste0(
                    "results/vsp_testing/", spf[k], "/model=vsp/predictions/", spf[k], "_model=vsp_mask_cog.tif"
                ))) next
                m <- rast(paste0(
                    "results/vsp_testing/", spf[k], "/model=vsp/predictions/", spf[k], "_model=vsp_mask_cog.tif"
                ))

                m <- m$fit_region_max_depth
                NAflag(m) <- 0

                r2 <- mask(r2, m)

                r1 <- mask(r1, r2)

                models_results [[mm]]  <- data.frame(
                    vsp_sub = spf[k],
                    Istat = dismo::nicheOverlap(
                    raster::raster(r1), raster::raster(r2), stat = "I"
                ),
                    model = models[mm]
                )
            }
            sub_results[[k]] <- dplyr::bind_rows(models_results)
        }
        results[[i]] <- dplyr::bind_rows(sub_results)
        results[[i]]$vsp <- species[i]
    }
    return(dplyr::bind_rows(results))
}

overlap <- check_overlap(c(1001, 1002, 1003, 1004))

overlap$vsp_sub <- gsub("_.$", "", overlap$vsp_sub)

overlap |> group_by(vsp, vsp_sub, model) |> summarise(mean = round(mean(Istat),1),sd=round(sd(Istat), 2), count = n()) |> 
    filter(model == "rf") |> filter(grepl("bias", vsp_sub))

##
library(terra)

sel_model = "rf"
r1 <- rast("data/virtual_species/key=1001/suitability_key1001.tif")
r2 <- rast("results/vsp_testing/taxonid=1001_occurrences_high_bias_1/model=vsp/predictions/taxonid=1001_occurrences_high_bias_1_model=vsp_method=rf_classification_ds_scen=current_cog.tif")

thresh <- arrow::read_parquet(
    "results/vsp_testing/taxonid=1001_occurrences_high_bias_1/model=vsp/metrics/taxonid=1001_occurrences_high_bias_1_model=vsp_what=thresholds.parquet"
)
thresh <- thresh |> filter(model == sel_model) |> select(p10) |> pull()

#r2[r2 < (thresh * 100)] <- 0
m <- rast("results/vsp_testing/taxonid=1001_occurrences_high_bias_1/model=vsp/predictions/taxonid=1001_occurrences_high_bias_1_model=vsp_mask_cog.tif")

m <- m$fit_region_max_depth
NAflag(m) <- 0

r2 <- mask(r2, m)

r1 <- mask(r1, r2)

plot(r2)

dismo::nicheOverlap(
    raster::raster(r1), raster::raster(r2), stat = "I"
)

