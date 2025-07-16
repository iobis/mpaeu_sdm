############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################### Main modelling function ############################


#' Model species distribution according to MPA Europe framework (species with
#' low number of records - ESM)
#'
#' @param species AphiaID of the species
#' @param group group of the species, according to the configuration file
#' @param species_dataset the path for the Parquet file containing species
#'   occurrences
#' @param outfolder the folder to output the results
#' @param outacro an unique acronym for the modelling
#' @param algorithms which algorithms to fit. For now only "esm" is available,
#'   and any other options will be ignored.
#' @param algo_opts options for the algorithms obtained through 
#'   [obissdm::sdm_options()]. IF `NULL` it will use the defaults. In this case,
#'   any attempt to set something with different 
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
model_species_esm <- function(species,
                              group,
                              species_dataset,
                              outfolder,
                              outacro,
                              algorithms = "esm",
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
  
  if (algorithms[1] != "esm") {
    cli::cli_alert_danger("For now only `algorithms` = 'esm' is available. Using {.val esm} instead of {.val {algorithms}}.")
    algorithms <- "esm"
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
  minptslim <- 10
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
    
    env <- .cm_check_coastal(species_data, env, coord_names, verb_2, filter_mode = "esm")
    
    # Restrict variables
    env$layers <- terra::subset(env$layers, env$hypothesis[[1]])
    env$hypothesis <- env$hypothesis[1]
    
    
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
    max_depth <- min(bath_pts)
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
      
      if (verb_2) cli::cli_alert_warning("Trying to get spatial blocks again with higher iteration value")

      sp_data <- mp_prepare_blocks(sp_data,
                                   method = "manual",
                                   block_types = "spatial_grid",
                                   n_iterate = 800,
                                   retry_if_zero = TRUE,
                                   manual_shp = block_grid,
                                   verbose = verb_2)
      
      if (any(table(sp_data$training$presence, sp_data$blocks$folds[["spatial_grid"]])[2,] == 0)) {
        stop("Blocks with less than 1 point. Failed.")
      }
    }
    vals_blocks <- table(sp_data$training$presence,
                         sp_data$blocks$folds[[1]])
    if(any(vals_blocks[2,] < 1)) {
      stop("Blocks with no presence point found. Failed.")
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
    model_log$model_details$hypothesis_tested <- env$hypothesis[1]
    model_log$model_details$best_hypothesis <- names(env$hypothesis)[1]
    model_log$model_details$variables <- env$hypothesis[[1]]
    
    # PART 5: FIT MODELS ----
    if (verb_1) cli::cli_alert_info("Starting model fitting")
    if (!is.null(algo_opts)) {
      if (any(!names(algo_opts) %in% names(obissdm::sdm_options("esm")))) {
        cli::cli_alert_info("Invalid `algo_opts`; using default.")
        algo_opts <- obissdm::sdm_options("esm")
      }
    } else {
      algo_opts <- obissdm::sdm_options("esm")
    }
    
    model_fits <- try(obissdm::sdm_module_esm(
      sp_data, options = algo_opts, verbose = verb_2
    ))
    
    treg <- obissdm::.get_time(treg, "Model fit")
    
    
    if (!inherits(model_fits, "try-error")) {
      if (verb_1) cli::cli_alert_info("Model fitting concluded.")
     
      # PART 6: PREDICTION ----
      if (length(model_fits) > 0) {
        avg_vars <- lapply(model_fits, function(x){
          cvm <- x$cv_metrics
          cvm <- cvm[[tg_metric]]
          if (sum(is.na(cvm)) >= 3) {
            0
          } else {
            mean(cvm, na.rm = T)
          }
        })
        avg_vars <- weighted.mean(unlist(avg_vars), attr(model_fits, "scores"))
      } else {
        avg_vars <- 0
      }
      
      if (avg_vars >= tg_threshold) {
        if (verb_1) cli::cli_alert_info("Good ESM, predicting.")
        # Predict models
        model_predictions <- lapply(seq_len(length(model_fits)), function(id) {
          
          model_name <- paste0("esm_part_", sprintf("%02d", c(seq_len(length(model_fits)))[id]))
          best_hyp <- names(env$hypothesis)
          
          if (verb_1) cli::cli_alert_info("Predicting model {id} - {model_name}")
          
          if (!dir.exists(pred_out)) fs::dir_create(pred_out)
          
          outmess <- paste0(pred_out, "taxonid=", species, "_model=", outacro, "_what=mess.tif")
          outshape <- paste0(pred_out, "taxonid=", species, "_model=", outacro, "_what=shape.tif")
          
          if (!file.exists(outmess)) {do_shape <- do_mess <- TRUE} else {do_shape <- do_mess <- FALSE} 
          
          scenarios <- data.frame(
            scenario = c("current", rep(c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
                                        each = 2)),
            year = c(NA, rep(c("dec50", "dec100"), 5))
          )
          # Reduced for testings
          # scenarios <- data.frame(scenario = c("current", "ssp126"), year = c(NA,"dec50"))
          # scenarios <- data.frame(scenario = c("current"), year = c(NA))
          
          for (k in seq_len(nrow(scenarios))) {
            
            if (is.na(scenarios$year[k])) {
              period <- NULL
            } else {
              period <- scenarios$year[k]
            }
            
            if (verb_1) cli::cli_alert_info("Predicting scenario {k} of {nrow(scenarios)}.")
            outpred <- paste0(pred_out, "taxonid=", species, "_model=", outacro,
                              "_method=", model_name, "_scen=", scenarios$scenario[k],
                              ifelse(is.null(period), "", paste0("_", period)), ".tif")
            
            env_to_pred <- obissdm::get_envofgroup(group,
                                                   depth = hab_depth, load_all = F,
                                                   scenario = scenarios$scenario[k],
                                                   period = period,
                                                   hypothesis = best_hyp,
                                                   env_folder = "data/env",
                                                   conf_file = "sdm_conf.yml", 
                                                   verbose = verb_2)
            env_to_pred <- terra::subset(env_to_pred, colnames(sp_data$training)[-1])

            if (best_hyp == "coastal") {
              env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
            }
            
            pred <- predict(model_fits[[id]], env_to_pred)
            
            names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))
            
            pred <- pred * 100
            pred <- as.int(pred)
            writeRaster(pred, outpred, overwrite = T, datatype = "INT1U")
            
            if (k == 1) {
              pred_f <- pred
            }
            
            if (do_mess) {
              # Save MESS
              if (verb_1) cli::cli_alert_info("Generating MESS map.")
              to_mess <- terra::aggregate(env_to_pred, 12)
              mess_map <- ecospat::ecospat.mess(
                na.omit(as.data.frame(to_mess, xy = T)),
                cbind(sp_data$coord_training, sp_data$training[,2:ncol(sp_data$training)]))
              mess_map_t <- to_mess[[1]]
              mess_map_t[] <- NA
              mess_map_t[cellFromXY(mess_map_t, mess_map[,1:2])] <- mess_map[,5]
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
              shape_data <- shape_data[c(which_p, which_a),]
              shape_data_coords <- sp_data$coord_training[c(which_p, which_a),]
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
                  outfile_fig <- paste0(outpath_fig,
                                        "taxonid=", species, "_model=", outacro, "_method=esm",
                                        "_what=shape.png")
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
                shape_res <- as.int(shape_res)
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
            pred_mess <- as.int(pred_mess)
            writeRaster(pred_mess, outmess, overwrite = T, datatype = "INT1U")
            cogeo_optim(outmess)
          }
          if (do_shape) {
            writeRaster(pred_shape, outshape, overwrite = T, datatype = "INT2U")
            cogeo_optim(outshape)
          }
          
          return(pred_f)
          
        })
        
        treg <- obissdm::.get_time(treg, "Model prediction")
        
        model_log$model_good <- 1
        model_log$model_good_metric <- tg_metric
        model_log$model_good_threshold <- tg_threshold
        
        
        
        # PART 7: GET RESPONSE CURVES ----
        if (verb_1) cli::cli_alert_info("Getting response curves")
        
        model_respcurves_all <- lapply(model_fits, function(model){
          vars <- model$variables
          sp_data_wv <- sp_data
          sp_data_wv$training <- sp_data_wv$training[,c("presence", vars)]
          obissdm::resp_curves(model, terra::subset(env$layers, vars), 
                               sp_data)
        })
        model_respcurves_all <- lapply(model_respcurves_all, function(x){
          x$value <- rep(1:100, (nrow(x)/100))
          x
        })
        model_respcurves_all <- do.call("rbind", model_respcurves_all)
        model_respcurves <- model_respcurves_all %>%
          group_by(variable, value) %>%
          summarise(response = mean(response),
                    base = mean(base), in_range = max(in_range)) %>%
          select(variable, response, base, in_range)
        class(model_respcurves) <- c("sdm_respcur", class(model_respcurves))
        
        if (!dir.exists(metric_out)) fs::dir_create(metric_out)
        model_name <- "esm"
        outfile <- paste0(metric_out,
                          "taxonid=", species, "_model=", outacro, "_method=", model_name,
                          "_what=respcurves.parquet")
        write_parquet(model_respcurves, outfile)
        
        if (!dir.exists(fig_out)) fs::dir_create(fig_out)
        outfile <- paste0(fig_out,
                          "taxonid=", species, "_model=", outacro, "_method=esm",
                          "_what=responsecurves.png")
        p <- plot(model_respcurves) + ggtitle(paste("AphiaID:", species, "| Model:", outacro, "| Algo: ESM"))
        ggplot2::ggsave(filename = outfile, plot = p,
                        width = 10, height = 8)
        
        
        treg <- obissdm::.get_time(treg, "Response curves")
        
        
        # Variables importance
        model_varimport <- obissdm::variable_importance_esm(
          model_fits, sp_data
        )
        outfile <- paste0(metric_out,
                          "taxonid=", species, "_model=", outacro, "_method=", model_name,
                          "_what=varimportance.parquet")
        
        arrow::write_parquet(model_varimport, outfile)
        
        treg <- obissdm::.get_time(treg, "Variable importance")
        
        # PART 8: ENSEMBLE ----
        to_ensemble <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/"),
                                  pattern = "\\.tif$", full.names = T)
        to_ensemble <- to_ensemble[grepl("_part_", to_ensemble)]
        
        to_ensemble_g <- stringr::str_extract(to_ensemble, "scen=*.*")
        to_ensemble_g <- gsub("scen=", "", to_ensemble_g)
        to_ensemble_g <- gsub("\\.tif", "", to_ensemble_g)
        to_ensemble_g <- unique(to_ensemble_g)
        
        for (un in seq_along(to_ensemble_g)) {
          to_ensemble_r <- rast(to_ensemble[grepl(to_ensemble_g[un], to_ensemble)])
          
          scores_ens <- attr(model_fits, "scores")
          
          ensemble_mean <- weighted.mean(to_ensemble_r, scores_ens, na.rm = T)
          
          if (grepl("current", to_ensemble_g[un])) {
            ensemble_curr <- ensemble_mean
          }
          
          names(ensemble_mean) <- gsub("_$", "", to_ensemble_g[un])
          
          outens <- paste0(pred_out, 
                           "taxonid=", species, "_model=", outacro, "_method=esm_scen=", gsub("_$", "", to_ensemble_g[un]), ".tif")
          
          ensemble_mean <- as.int(ensemble_mean)
          writeRaster(ensemble_mean, outens,
                      overwrite = T, datatype = "INT1U")
          res_ens <- cogeo_optim(outens)
          rm(ensemble_mean)
        }
        
        fs::file_delete(to_ensemble)
        
        model_log$model_result$esm <- "succeeded"
        
        # Evaluate ensemble
        ensemble_curr <- ensemble_curr / 100
        
        fit_ens <- terra::extract(ensemble_curr, sp_data$coord_training, ID = F)
        
        fit_ens_eval <- obissdm::eval_metrics(sp_data$training$presence, fit_ens[,1])
        
        if (!is.null(sp_data$eval_data)) {
          eval_ens <- terra::extract(ensemble_curr, sp_data$coord_eval, ID = F)
          
          eval_ens_eval <- obissdm::eval_metrics(sp_data$eval_data$presence, eval_ens[,1])
          
          ensemble_metrics <- rbind(cbind(as.data.frame(t(fit_ens_eval)), what = "full", origin = "fit_data"),
                                    cbind(as.data.frame(t(eval_ens_eval)), what = "full", origin = "eval_data"))
        } else {
          ensemble_metrics <- cbind(as.data.frame(t(fit_ens_eval)), what = "full", origin = "fit_data")
        }
        
        # Get list of metrics
        metrics_list <- lapply(model_fits, function(imod){
          imod$cv_metrics
        })
        
        # Convert the list of data frames into a 3-dimensional array
        scores <- attr(model_fits, "scores")
        metrics_array <- array(unlist(metrics_list), dim = c(nrow(metrics_list[[1]]), ncol(metrics_list[[1]]), length(metrics_list)))
        
        # Calculate the mean across the third dimension (which represents the data frames)
        metrics_mean_df <- apply(metrics_array, c(1, 2), weighted.mean, w = scores, na.rm = T)
        metrics_mean_sd <- apply(metrics_array, c(1, 2), sd, na.rm = T)
        
        # Convert the resulting matrix back to a data frame and set the column names
        metrics_mean_df <- as.data.frame(metrics_mean_df)
        colnames(metrics_mean_df) <- colnames(metrics_list[[1]])
        metrics_mean_df$what <- "mean"
        metrics_mean_df$origin <- "avg_fit"
        
        metrics_mean_sd <- as.data.frame(metrics_mean_sd)
        colnames(metrics_mean_sd) <- colnames(metrics_list[[1]])
        metrics_mean_sd$what <- "sd"
        metrics_mean_sd$origin <- "avg_fit"
        
        avg_fit_metrics <- rbind(metrics_mean_df, metrics_mean_sd)
        
        if (!is.null(sp_data$eval_data)) {
          
          # Get list of metrics
          metrics_list <- lapply(good_models, function(id){
            model_fits[[id]]$eval_metrics
          })
          metrics_list <- do.call("rbind", metrics_list)
          
          # Calculate the mean across
          eval_metrics_mean_df <- apply(metrics_list, 2, weighted.mean, w = scores, na.rm = T)
          eval_metrics_mean_sd <- apply(metrics_list, 2, sd, na.rm = T)
          
          # Convert the resulting matrix back to a data frame
          eval_metrics_mean_df <- as.data.frame(t(eval_metrics_mean_df))
          eval_metrics_mean_df$what <- "mean"
          eval_metrics_mean_df$origin <- "avg_eval"
          
          eval_metrics_mean_sd <- as.data.frame(t(eval_metrics_mean_sd))
          eval_metrics_mean_sd$what <- "sd"
          eval_metrics_mean_sd$origin <- "avg_eval"
          
          avg_fit_metrics <- rbind(avg_fit_metrics, eval_metrics_mean_df, eval_metrics_mean_sd)
          
        }
        
        ensemble_metrics <- rbind(ensemble_metrics, avg_fit_metrics)
        
        arrow::write_parquet(ensemble_metrics,
                             paste0(metric_out, 
                                    "taxonid=", species, "_model=", outacro,
                                    "_method=esm_what=cvmetrics.parquet"))
        
        # PART 9: SAVE BINARIZATION INFO ----
        if (verb_1) cli::cli_alert_info("Saving binarization info")
        # P10 threshold
        thresh_p10_mtp <- lapply(ensemble_curr, function(pred) {
          predv <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1,], ID = F)[,1]
          p10 <- ceiling(length(predv) * 0.9)
          data.frame(p10 = rev(sort(predv))[p10], mtp = min(na.omit(predv)))
        })
        
        # Max Specificity + Sensitivity
        thresh_maxss <- lapply(ensemble_curr, function(pred) {
          pred_p <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1,], ID = F)[,1]
          pred_a <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 0,], ID = F)[,1]
          
          pred_eval <- predicts::pa_evaluate(p = pred_p, a = pred_a)
          
          pred_thresh <- predicts::threshold(pred_eval)
          
          pred_thresh
        })
        
        thresh_p10_mtp <- bind_rows(thresh_p10_mtp)
        thresh_maxss <- bind_rows(thresh_maxss)
        
        thresh <- cbind(thresh_p10_mtp, thresh_maxss)
        thresh$model <- "esm"
        
        arrow::write_parquet(thresh, paste0(metric_out,
                                     "taxonid=", species, "_model=", outacro,
                                     "_what=thresholds.parquet"))
        
        
        
        # PART 10: CREATE MASKS AND SAVE ----
        if (verb_1) cli::cli_alert_info("Saving masks")
        best_hyp <- names(env$hypothesis)
        .cm_save_masks(
          ecoreg_occ, ecoreg_sel, multi_mod = list(best_model = best_hyp),
          env, model_predictions, sp_data,
          pred_out, species, outacro, coord_names, max_depth
        )
        treg <- obissdm::.get_time(treg, "Masks")
        
        
        
        # PART 11: EVALUATE MODELS USING OTHER TECHNIQUES ----
        if (verb_1) cli::cli_alert_info("Performing post-evaluation")
        model_predictions <- list(ensemble_curr)
        good_models <- 1
        
        model_log <- .cm_posteval_all(
          post_eval, model_predictions, sp_data,
          thresh_p10_mtp, algorithms, hab_depth,
          good_models, model_log, metric_out,
          species, outacro,
          list(model_varimport), env, verb_1
        )
        treg <- obissdm::.get_time(treg, "Post-evaluation")
        
        
        # PART 12: SAVE MODELS OBJECTS ----
        if (verb_1) cli::cli_alert_info("Saving models objects")
        
        # Update log with best parameters
        model_log$model_bestparams$esm <- list(
           individual_parameters = bind_rows(
             lapply(
               seq_len(length(model_fits)), function(x) cbind(model_fits[[x]]$parameters, part = x)
             )
           ),
           scores = attr(model_fits, "scores"),
           variable_combinations = lapply(
             seq_len(length(model_fits)), function(x) model_fits[[x]]$variables
           )
        )
        
        # Save models
        model_fits <- lapply(model_fits, function(x){
          reduced_model <- try(.reduce_model_file(x))
          if (!inherits(reduced_model, "try-error")) {
            reduced_model <- x
          }
          reduced_model
        })
        class(model_fits) <- c("sdm_esm_result", class(model_fits))
        
        outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
        if (!dir.exists(outpath)) fs::dir_create(outpath)
        outfile <- paste0(outpath,
                          "taxonid=", species, "_model=", outacro, "_method=esm",
                          "_what=model.rds")
        
        saveRDS(model_fits, file = outfile)
        
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
        if (verb_1) cli::cli_alert_warning("Average ESM metrics bad")
        st_status <- "no_good_model"
      }
      
    } else {
      if (verb_1) cli::cli_alert_warning("ESM failed")
      st_status <- "all_failed"
      if (verb_2) print(model_fits)
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
