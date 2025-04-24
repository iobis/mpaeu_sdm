############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################### Main modelling function ############################


#' Model species distribution according to MPA Europe framework
#'
#' @param species AphiaID of the species
#' @param group group of the species, according to the configuration file
#' @param species_dataset the path for the Parquet file containing species
#'   occurrences
#' @param outfolder the folder to output the results
#' @param outacro an unique acronym for the modelling
#' @param algorithms which algorithms to fit. See [obissdm::sdm_options()] for a
#'   list of available algorithms
#' @param algo_opts options for the algorithms obtained through 
#'   [obissdm::sdm_options()]. IF `NULL` it will use the defaults
#' @param limit_by_depth should the study area be limited by depth? If so, the
#'   depth is extracted from occurrence records and after applying the buffer
#'   the layers are masked. Recommended to keep `TRUE`
#' @param depth_buffer the depth buffer to be added in case `limit_by_depth =
#'   TRUE`. Can be 0 to apply no buffer (not recommended)
#' @param assess_bias if `TRUE`, perform tests for potential spatial bias on
#'   occurrence records
#' @param post_eval character vector with names of post-evaluation methods to 
#'   try. Should be one of `sst` (thermal niche), `niche` (niche equivalency)
#'   or `hyper` (comparison with hypervolume niche estimation)
#' @param tg_metric target metric to be used to assess model quality
#' @param tg_threshold threshold for the target metric to a model be considered
#'   of good quality
#' @param quad_samp number of quadrature points to be sampled. Can be also a
#'   number between 0 (higher than) and 1, expressing the percentage of available
#'   environmental cells to be sampled as quadrature points
#' @param cleanup if `TRUE` (recommended), if the folder for the species already
#'   exists, it will be removed
#' @param verbose if `TRUE` print essential messages. Can also be numeric: 0 for
#'   no messages, 1 for progress messages and 2 for all messages.
#'
#' @return nothing, saved files
#' @export
#' 
#' @details
#' This function is a wrapper around functions from the [obissdm] package, to
#' provide a consistent pipeline for generating the SDMs according to the 
#' MPA Europe project.
#' 
#' Output will be in the format:
#' {outfolder}/taxonid={AphiaID}/model={outacro}/predictions OR models OR metrics/
#' 
#' - predictions: contains all raster files, optimized for cloud (COG)
#' - models: contain model object in zip files
#' - metrics: contain relevant metrics
#' 
#' A single `.json` file is saved on the root of the model containing details
#' of the model fit. Also a `.parquet` file is saved with the occurrence data
#' used to fit the model (i.e. already filtered and standardized)
#'
#' @examples
#' \dontrun{
#' model_species(124316, "benthic_invertebrates", "sp_dataset.parquet",
#' "results", "mr1")
#' }
model_species <- function(species,
                          group,
                          species_dataset,
                          outfolder,
                          outacro,
                          algorithms = c("maxent", "rf", "brt", "lasso"),
                          algo_opts = NULL,
                          limit_by_depth = TRUE,
                          depth_buffer = 500,
                          assess_bias = TRUE,
                          post_eval = c("sst", "niche", "hyper"),
                          tg_metric = "cbi",
                          tg_threshold = 0.3,
                          quad_samp = 50000,
                          cleanup = TRUE,
                          max_mem = 0.6,
                          verbose = FALSE) {
  
  terra::terraOptions(memfrac = max_mem)

  # Check verbosity
  verb_1 <- verb_2 <- FALSE
  if (is.numeric(verbose)) {
    if (verbose == 1) {
      verb_1 <- TRUE
    } else if (verbose == 2) {
      verb_1 <- verb_2 <- TRUE
    }
  } else {
    if (verbose) {
      verb_1 <- TRUE
    }
  }
  
  if (verb_1) cli::cli_alert_info("Starting model for species {species}")

  if (cleanup && dir.exists(file.path(outfolder, paste0("taxonid=", species)))) {
    fs::dir_delete(file.path(outfolder, paste0("taxonid=", species)))
  }
  
  # Record timing
  treg <- obissdm::.get_time()
  
  # Define global objects and fill log
  coord_names <- c("decimalLongitude", "decimalLatitude")
  model_log <- obissdm::gen_log(algorithms)
  model_log$taxonID <- species
  model_log$group <- group
  model_log$model_acro <- outacro
  model_log$model_date <- Sys.Date()
  if (!is.null(algo_opts)) {
    algo_opts_nams <- names(algo_opts)
    for (opt in algo_opts_nams) model_log$algorithms_parameters[[opt]] <- algo_opts[[opt]]
  }
  
  # Define output paths
  metric_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
  pred_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/")
  mod_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
  fig_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/figures/")
  
  # Load species data
  if (verb_1) cli::cli_alert_info("Reading data")
  species_data <- .cm_load_species(species, species_dataset)
  
  treg <- obissdm::.get_time(treg, "Species data loading")
  
  # If data is available, proceeds
  minptslim <- 30
  if (nrow(species_data) > 0 && nrow(species_data[species_data$data_type == "fit_points",]) >= minptslim) {
    
    if (verb_1) cli::cli_alert_info("Enough number of points, proceeding")
    # PART 1: DATA LOADING ----
    # Load species information
    species_name <- species_data$species[1]
    model_log$scientificName <- species_name
    
    # Load and select ecological information for the species
    hab <- .cm_check_habitat(species)
    hab_depth <- hab_to_depth(hab)
    model_log$hab_depth <- hab_depth
    
    # Load environmental data
    env <- get_envofgroup(group, depth = hab_depth, load_all = T,
                          conf_file = "sdm_conf.yml", verbose = verb_2)
    
    env <- .cm_check_coastal(species_data, env, coord_names, verb_2)
    
    
    # Load Realms shapefile (we call ecoregions for convenience)
    ecoregions <- vect("data/shapefiles/MarineRealms_BO.shp")
    
    # Split the dataset
    fit_pts <- split_ds(species_data)
    eval_pts <- split_ds(species_data, "eval")
    model_log$n_init_points <- nrow(fit_pts)
    
    treg <- obissdm::.get_time(treg, "Data loading")
    
    
    
    # PART 2: DATA PREPARING ----
    if (verb_1) cli::cli_alert_info("Preparing data")

    # Check which eco-regions are covered by the points
    prep_eco <- .cm_check_ecoregions(
      ecoregions, fit_pts, eval_pts, env, limit_by_depth, depth_buffer,
      coord_names, model_log, verb_1
    )
    env <- prep_eco$env
    model_log <- prep_eco$model_log
    bath_pts <- prep_eco$bath_pts
    ecoreg_occ <- prep_eco$ecoreg_occ
    ecoreg_sel <- prep_eco$ecoreg_sel
    rm(prep_eco)
    
    # Check if species is from shallow/coastal areas and remove distcoast/bathymetry
    env <- .cm_check_shallow(bath_pts, env, verb_1)
    
    # Assess SAC
    fit_pts_sac <- try(obissdm::outqc_sac_mantel(fit_pts, 
                                                 env_layers = terra::subset(env$layers, env$hypothesis[[1]]),
                                                 plot_result = FALSE,
                                                 verbose = verb_2))
    
    if (inherits(fit_pts_sac, "try-error")) {
      fit_pts_sac <- fit_pts
    } 

    # Make data object
    quad_n <- .cm_calc_quad(env, quad_samp, fit_pts_sac)
    model_log$model_details$background_size <- quad_n 

    sp_data <- mp_prepare_data(fit_pts_sac, eval_data = eval_pts,
                               species_id = species_name,
                               env_layers = env$layers,
                               quad_number = quad_n,
                               verbose = verb_2)

    block_grid <- get_block_grid(sp_data, env$layers,
                                 sel_vars = env$hypothesis$basevars,
                                 verbose = verb_2)

    sp_data <- mp_prepare_blocks(sp_data,
                                 method = "manual",
                                 block_types = "spatial_grid",
                                 n_iterate = 300,
                                 retry_if_zero = TRUE,
                                 manual_shp = block_grid,
                                 verbose = verb_2)
    
    if (any(table(sp_data$training$presence, sp_data$blocks$folds[["spatial_grid"]])[2,] == 0)) {
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
      multi_mod <- sdm_multhypo(sdm_data = sp_data, sdm_method = "rf",
                                    variables = env$hypothesis,
                                    options = algo_opts[["rf"]],
                                    verbose = verb_2)
    } else {
      if (verb_1) cli::cli_alert_info("Only one hypothesis available - skipping")
      multi_mod <- list(best_model = names(env$hypothesis),
                            best_variables = env$hypothesis[[1]])
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
    fit_obj <- .cm_model_fit(algorithms, algo_opts, sp_data,
                             model_log, verb_1, verb_2)
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
        model_predictions <- .cm_predict_models(good_models, model_fits, 
                                                multi_mod, pred_out, group, 
                                                hab_depth, sp_data, outfolder, 
                                                outacro, species, env,
                                                verb_1, verb_2)
        
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
        thresh_p10_mtp <- .cm_save_bin_info(model_predictions, sp_data, 
                                            algorithms, good_models,
                                            metric_out, species, outacro)


        # PART 10: CREATE MASKS AND SAVE ----
        if (verb_1) cli::cli_alert_info("Saving masks")
        .cm_save_masks(
          ecoreg_occ, ecoreg_sel, multi_mod,
          env, model_predictions, sp_data,
          pred_out, species, outacro, coord_names
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
        .cm_save_metrics(model_fits, metric_out, species,
                         outacro, good_models)
        
        # Update log with best parameters
        for (al in 1:length(algorithms)) {
          if (!is.null(model_fits[[al]])) {
            model_log$model_bestparams[[al]] <- model_fits[[al]]$parameters
          }
        }
        
        # Save models
        .cm_save_models(model_fits, outfolder, species, outacro, good_models)
        
        # Save fit points
        arrow::write_parquet(sp_data$coord_training[sp_data$training$presence == 1,],
                      paste0(outfolder, "/taxonid=", species, "/model=", outacro,
                             "/taxonid=", species, "_model=", outacro, "_",
                             "what=fitocc.parquet"))
        
        treg <- as.data.frame(unclass(treg))
        treg <- data.frame(what = row.names(treg), time_mins = treg[,1])
        model_log$timings <- treg
        
        outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/")
        outfile <- paste0(outpath,
                          "taxonid=", species, "_model=", outacro, "_what=log.json")
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
    
  } else {
    if (nrow(species_data) > 0 && nrow(species_data) < minptslim) {
      if (verb_1) cli::cli_alert_warning("Low number of points, failed")
      st_status <- "low_data"
    } else {
      if (verb_1) cli::cli_alert_warning("No data available, failed")
      st_status <- "no_data"
    }
  }
  
  return(st_status)
  
}
