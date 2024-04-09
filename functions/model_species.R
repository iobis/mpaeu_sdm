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
#' @param limit_by_depth should the study area be limited by depth? If so, the
#'   depth is extracted from occurrence records and after applying the buffer
#'   the layers are masked. Recommended to keep `TRUE`
#' @param depth_buffer the depth buffer to be added in case `limit_by_depth =
#'   TRUE`. Can be 0 to apply no buffer (not recommended)
#' @param assess_bias if `TRUE`, perform tests for potential spatial bias on
#'   occurrence records
#' @param correct_bias if `TRUE`, apply a bias layer to try to correct the
#'   spatial bias
#' @param verbose if `TRUE` print messages
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
#' of the model fit.
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
                          limit_by_depth = TRUE,
                          depth_buffer = 500,
                          assess_bias = TRUE,
                          correct_bias = TRUE,
                          verbose = FALSE) {
  
  if (verbose) cli::cli_alert_info("Starting model for species {species}")
  
  # Record timing
  treg <- obissdm::.get_time()
  
  # Define global objects and fill log
  coord_names <- c("decimalLongitude", "decimalLatitude")
  model_log <- obissdm::gen_log(algorithms)
  model_log$taxonID <- species
  model_log$model_acro <- outacro
  model_log$model_date <- Sys.Date()
  
  # Load species data
  if (verbose) cli::cli_alert_info("Reading data")
  if (tools::file_ext(species_dataset) == "parquet") {
    if (utils::file_test("-d", species_dataset)) {
      ds <- pl$scan_parquet(file.path(species_dataset, "**/*"))
    } else {
      ds <- pl$scan_parquet(species_dataset)
    }
    species_data <- ds$filter(pl$col("taxonID") == species)$collect()
    species_data <- species_data$to_data_frame()
  } else if (tools::file_ext(species_dataset) %in% c("sqlite", "db")) {
    #to do!
  } else if (any(grepl("key", list.files(species_dataset, recursive = T)))) {
    lf <- list.files(species_dataset, recursive = T)
    lf <- lf[grepl(paste0("key=", species), lf)]
    if (length(lf) > 0) {
      species_data <- arrow::read_parquet(paste0(species_dataset, "/", lf))
    } else {
      species_data <- matrix(nrow = 0, ncol = 1)
    }
  } else {
    cli::cli_abort("No folder or files were found in {.path {species_dataset}}")
  }
  
  treg <- obissdm::.get_time(treg, "Species data loading")
  
  # If data is available, proceeds
  if (nrow(species_data) >= 15) {
    
    if (verbose) cli::cli_alert_info("Enough number of points, proceeding")
    # PART 1: DATA LOADING ----
    # Load species information
    species_name <- species_data$species[1]
    model_log$scientificName <- species_name
    
    # Load ecological information
    eco_info <- open_csv_dataset("data/species_ecoinfo.csv") %>%
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
    
    # Load environmental data
    env <- get_envofgroup(group, depth = hab_depth, load_all = T,
                          conf_file = "sdm_conf.yml", verbose = verbose)
    
    if ("coastal" %in% names(env$hypothesis)) {
      # Test if is within europe/coastal
      europe_starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
      europe_starea <- ext(europe_starea)
      
      is_within <- is.related(vect(species_data[,coord_names], geom = coord_names),
                              europe_starea, "intersects")
      if (!all(is_within)) {
        env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
        env$layers <- subset(env$layers, "wavefetch", negate = T)
        if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
      } else {
        
        wave <- rast("data/env/terrain/wavefetch.tif")
        wave_pts <- extract(wave, species_data[,coord_names])
        
        if (sum(is.na(wave_pts[,2])) > ceiling(nrow(species_data) * .05)) {
          env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
          if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
        } else if ((nrow(species_data) - sum(is.na(wave_pts[,2]))) < 15) {
          env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
          if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
        } else {
          env$layers <- crop(env$layers, europe_starea)
        }
      }
    }
    
    
    # Load ecoregions shapefile
    ecoregions <- vect("data/shapefiles/MarineRealms.shp")
    
    # Split the dataset
    fit_pts <- split_ds(species_data)
    eval_pts <- split_ds(species_data, "eval")
    
    treg <- obissdm::.get_time(treg, "Data loading")
    
    
    
    # PART 2: DATA PREPARING ----
    if (verbose) cli::cli_alert_info("Preparing data")
    # Check which eco-regions are covered by the points
    ecoreg_unique <- ecoregions$Realm[is.related(ecoregions,
                                                 vect(bind_rows(fit_pts, eval_pts), 
                                                      geom = coord_names), 
                                                 "intersects")]
    
    model_log$model_details$ecoregions <- ecoreg_unique
    ecoreg_occ <- ecoregions[ecoregions$Realm %in% ecoreg_unique,]
    
    # Apply a buffer to ensure that all areas are covered
    sf::sf_use_s2(FALSE)
    ecoreg_occ_buff <- suppressMessages(
      suppressWarnings(vect(sf::st_buffer(sf::st_as_sf(ecoreg_occ), 0.2)))
    )
    
    adj_ecoreg <- ecoregions$Realm[is.related(ecoregions, ecoreg_occ_buff,
                                              "intersects")]
    
    # Mask areas
    ecoreg_sel <- ecoregions[ecoregions$Realm %in% unique(
      c(ecoreg_unique, adj_ecoreg)
    ),]
    model_log$model_details$ecoregions_included <- unique(ecoreg_sel$Realm)
    
    # Apply a buffer to ensure that all areas are covered
    sf::sf_use_s2(FALSE)
    ecoreg_sel <- suppressMessages(
      suppressWarnings(vect(sf::st_buffer(sf::st_as_sf(ecoreg_sel), 0.5)))
    )
    
    # Limit to max extension based on distance to the points
    max_ext <- ext(vect(rbind(fit_pts, eval_pts), geom = coord_names))
    max_ext <- ext(as.vector(max_ext) + c(-20, 20, -20, 20))
    ecoreg_sel <- crop(ecoreg_sel, max_ext)
    
    # Limit by depth if TRUE
    if (limit_by_depth) {
      if (verbose) cli::cli_alert_info("Limiting by depth")
      # To decide if this step will be kept:
      # Load bathymetry layer
      bath <- rast("data/env/terrain/bathymetry_mean.tif")
      bath <- mask(crop(bath, ecoreg_sel), ecoreg_sel)
      
      bath_pts <- terra::extract(bath, bind_rows(fit_pts, eval_pts))
      
      bath_range <- range(bath_pts[,2])
      bath_range[1] <- bath_range[1] - depth_buffer
      bath_range[2] <- ifelse((bath_range[2] + depth_buffer) > 0, 
                              0, bath_range[2] + depth_buffer)
      
      bath[bath < bath_range[1] | bath > bath_range[2]] <- NA
      
      if ("coastal" %in% names(env$hypothesis)) {
        bath <- crop(bath, europe_starea)
        env$layers <- mask(crop(env$layers, ecoreg_sel), bath)
        env$layers <- mask(env$layers, env$layers$wavefetch)
      } else {
        env$layers <- mask(crop(env$layers, ecoreg_sel), bath)
      }
      
      model_log$model_details$limited_by_depth <- TRUE
      model_log$model_details$depth_buffer <- depth_buffer
      
    } else {
      if ("coastal" %in% names(env$hypothesis)) {
        ecoreg_sel <- crop(ecoreg_sel, europe_starea)
        env$layers <- mask(env$layers, ecoreg_sel)
      } else {
        env$layers <- mask(env$layers, ecoreg_sel)
      }
    }
    
    # Assess SAC
    fit_pts_sac <- try(obissdm::outqc_sac(fit_pts, 
                                          env_layers = subset(env$layers, env$hypothesis[[1]]),
                                          plot_result = FALSE,
                                          verbose = verbose))
    
    if (inherits(fit_pts_sac, "try-error")) {
      fit_pts_sac <- fit_pts
    } 
    
    # Make data object
    quad_n <- ifelse(nrow(as.data.frame(env$layers, xy = T, na.rm = T)) < 50000,
                     round(nrow(as.data.frame(env$layers, xy = T, na.rm = T))), 50000)
    
    sp_data <- mp_prepare_data(fit_pts_sac, eval_data = eval_pts,
                               species_id = species_name,
                               env_layers = env$layers,
                               quad_number = quad_n,
                               verbose = verbose)
    
    block_grid <- get_block_grid(sp_data, env$layers,
                                 sel_vars = env$hypothesis$basevars,
                                 verbose = verbose)
    
    sp_data <- mp_prepare_blocks(sp_data,
                                 method = "manual",
                                 block_types = "spatial_grid",
                                 n_iterate = 300,
                                 retry_if_zero = TRUE,
                                 manual_shp = block_grid,
                                 verbose = verbose)
    
    model_log$model_details$block_size <- unname(sp_data$blocks$grid_resolution)
    model_log$model_fit_points <- sum(sp_data$training$presence)
    model_log$model_eval_points <- sum(sp_data$eval_data$presence)
    
    treg <- obissdm::.get_time(treg, "Data preparing")
    
    # PART 3: ASSESS SPATIAL BIAS ----
    if (assess_bias) {
      if (verbose) cli::cli_alert_info("Assessing spatial bias")
      # Add some way of assessing spatial bias
      model_log$model_details$control_bias <- correct_bias
      
      # Get spatstat statistics and point process for control
      spat_im <- as.im(as.data.frame(aggregate(env$layers[[1]], 10, na.rm = T), xy = T))
      spat_window <- as.owin(spat_im)
      # Define ppp object based on  point locations
      spat_ppp <- ppp(x = sp_data$coord_training$decimalLongitude[sp_data$training$presence == 1],
                      y = sp_data$coord_training$decimalLatitude[sp_data$training$presence == 1],
                      window = spat_window)
      # calculate envelope around L-hat estimates.
      spat_env_k <- envelope(spat_ppp, Kest, verbose = F)
      spat_env_l <- envelope(spat_ppp, Lest, verbose = F)
      
      # Point process model (?)
      spat_quad <- ppp(x = sp_data$coord_training$decimalLongitude[sp_data$training$presence == 0],
                       y = sp_data$coord_training$decimalLatitude[sp_data$training$presence == 0],
                       window = spat_window)
      spat_quad_sc <- quadscheme(data = spat_ppp, dummy = spat_quad, method = "dirichlet")
      spat_data_obj <- lapply(1:nlyr(env$layers), function(x){
        as.im(as.data.frame(aggregate(env$layers[[x]], 4, na.rm = T), xy = T))
      })
      names(spat_data_obj) <- names(env$layers)
      spat_ppm <- ppm(spat_quad_sc,
                      trend = as.formula(paste("~",
                                               paste(names(spat_data_obj)[
                                                 names(spat_data_obj) %in% env$hypothesis[[1]]
                                               ], collapse = "+"))),
                      covariates = spat_data_obj)
      
      spat_pred <- predict(spat_ppm, covariates = spat_data_obj, ngrid = c(spat_data_obj[[1]]$dim[1], spat_data_obj[[1]]$dim[2]))
      
      if (correct_bias) {
        bias_layer <- density(spat_ppp, sigma = 1)
        bias_layer <- rast(bias_layer)
        bias_layer <- disagg(bias_layer, 10)
        crs(bias_layer) <- crs("EPSG:4326")
        bias_layer <- project(bias_layer, env$layers[[1]])
        bias_layer <- mask(bias_layer, env$layers[[1]])
      }
        
      treg <- obissdm::.get_time(treg, "Sample bias assessment")
    }
    
    
    # PART 4: MODEL SELECTION ----
    if (verbose) cli::cli_alert_info("Model selection")
    multi_mod_max <- sdm_multhypo(sdm_data = sp_data, sdm_method = "maxent",
                                  variables = env$hypothesis,
                                  verbose = verbose)
    
    if (correct_bias) {
      if (verbose) cli::cli_alert_info("Testing bias correction")
      temp_data <- sp_data
      temp_data$training <- temp_data$training[, c("presence", multi_mod_max$best_variables)]
      temp_data$eval_data <- temp_data$eval_data[, c("presence", multi_mod_max$best_variables)]
      
      temp_data$training$biasgrid <- extract(bias_layer, temp_data$coord_training, ID = F)[,1]
      temp_data$eval_data$biasgrid <- extract(bias_layer, temp_data$coord_eval, ID = F)[,1]
      
      not_na <- !is.na(temp_data$training$biasgrid)
      not_na_eval <- !is.na(temp_data$eval_data$biasgrid)
      
      temp_data$training <- temp_data$training[not_na,]
      temp_data$coord_training <- temp_data$coord_training[not_na,]
      
      temp_data$eval_data <- temp_data$eval_data[not_na_eval,]
      temp_data$coord_eval <- temp_data$coord_eval[not_na_eval,]
      
      max_bias <- sdm_fit(temp_data)
      
      bias_metrics <- apply(max_bias$cv_metrics, 2, mean, na.rm = T)
      model_metrics <- apply(multi_mod_max$model$cv_metrics, 2, mean, na.rm = T)
      
      if (bias_metrics[["cbi"]] > (model_metrics[["cbi"]] + 0.1)) {
        # TODO: continue here
        # TODO: bias layer needs to be other, not from the single species
        # group or something - to see
      }
    }
    
    # Prepare data for the best model
    sp_data$training <- sp_data$training[, c("presence", multi_mod_max$best_variables)]
    sp_data$eval_data <- sp_data$eval_data[, c("presence", multi_mod_max$best_variables)]
    
    treg <- obissdm::.get_time(treg, "Model selection")
    
    
    # PART 4: FIT MODELS ----
    if (verbose) cli::cli_alert_info("Starting model fitting")
    model_fits <- list()
    for (i in 1:length(algorithms)) {
      if (verbose) cli::cli_alert_info("Fitting {algorithms[i]}")
      model_fits[[i]] <- try(do.call(paste0("sdm_module_", algorithms[i]), list(sp_data, verbose = verbose)))
      if (inherits(model_fits[[i]], "try-error")) {
        model_log$model_result[[algorithms[i]]] <- "failed"
        model_fits[[i]] <- NULL
      } else {
        model_log$model_result[[algorithms[i]]] <- "succeeded"
      }
    }
    
    # For some reason was not registering in the log... using for instead
    # model_fits <- lapply(algorithms, function(algo) {
    #   fit <- try(do.call(paste0("sdm_module_", algo), list(sp_data)))
    #   if (inherits(fit, "try-error")) {
    #     model_log$model_result[[algo]] <- "failed"
    #     return(NULL)
    #   } else {
    #     model_log$model_result[[algo]] <- "succeeded"
    #     return(fit)
    #   }
    # })
    
    treg <- obissdm::.get_time(treg, "Model fit")
    
    
    if (any(unlist(model_log$model_result) == "succeeded")) {
      if (verbose) cli::cli_alert_info("Model fitting concluded. Checking if models are good")
      # PART 5: PREDICTION ----
      
      # Define thresholds
      tg_metrics <- "cbi"
      tg_threshold <- 0.3
      
      # Get what models are above threshold
      good_models <- lapply(model_fits, function(model){
        if (!is.null(model)) {
          cv_res <- model$cv_metrics
          cv_res <- apply(cv_res, 2, mean, na.rm = T)
          the_metric <- cv_res[[tg_metrics]]
          if (!is.na(the_metric)) {
            if (the_metric >= tg_threshold) {
              if (sum(is.na(model$cv_metrics[[tg_metrics]])) >= 4) {
                return(FALSE)
              } else {
                return(TRUE)
              }
            } else {
              return(FALSE)
            }
          } else {
            return(FALSE)
          }
        } else {
          return(FALSE)
        }
      })
      
      good_models <- which(unlist(good_models))
      
      if (length(good_models) > 0) {
        if (verbose) cli::cli_alert_info("Good models available, predicting.")
        # Predict models
        model_predictions <- lapply(good_models, function(id) {
          model_name <- model_fits[[id]]$name
          
          best_hyp <- multi_mod_max$best_model
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/")
          
          if (!dir.exists(outpath)) fs::dir_create(outpath)
          
          outmess <- paste0(outpath, "taxonid=", species, "_model=", outacro, "_what=mess.tif")
          
          if (!file.exists(gsub("mess", "mess_cog", outmess))) {do_mess <- TRUE} else {do_mess <- FALSE} 
          
          scenarios <- data.frame(
            scenario = c("current", rep(c("ssp116", "ssp245", "ssp370", "ssp460", "ssp585"),
                                        each = 2)),
            year = c(NA, rep(c("dec50", "dec100"), 5))
          )
          
          for (k in 1:nrow(scenarios)) {
            
            if (is.na(scenarios$year[k])) {
              period <- NULL
            } else {
              period <- scenarios$year[k]
            }
            
            if (verbose) cli::cli_alert_info("Predicting scenario {k} of {nrow(scenarios)}.")
            
            env_to_pred <- obissdm::get_envofgroup(group,
                                                   depth = hab_depth, load_all = F,
                                                   scenario = scenarios$scenario[k],
                                                   period = period,
                                                   hypothesis = best_hyp,
                                                   env_folder = "data/env",
                                                   conf_file = "sdm_conf.yml", 
                                                   verbose = verbose)
            
            pred <- predict(model_fits[[id]], env_to_pred)
            
            names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))
            
            if (do_mess) {
              # Save MESS
              if (!exists("to_mess")) {
                maxpt <- nrow(as.data.frame(env_to_pred, xy = T, na.rm = T))
                to_mess <- sample(1:maxpt,
                                  ifelse(maxpt > 10000, 10000, maxpt))
              }
              mess_map <- ecospat::ecospat.mess(
                na.omit(as.data.frame(env_to_pred, xy = T)[to_mess,]), 
                cbind(sp_data$coord_training, sp_data$training[,2:ncol(sp_data$training)]))
              mess_map_t <- env_to_pred[[1]]
              mess_map_t[] <- NA
              mess_map_t[cellFromXY(mess_map_t, mess_map[,1:2])] <- mess_map[,5]
              mess_map <- mask(mess_map_t, env_to_pred[[1]])
              
              names(mess_map) <- names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))
              
              if (k == 1) {
                pred_mess <- mess_map
              } else {
                pred_mess <- c(pred_mess, mess_map)
              }
            }
            
            
            if (k == 1) {
              pred_f <- pred
            } else {
              pred_f <- c(pred_f, pred)
            }
          }
          
          writeRaster(pred_f, paste0(outpath, "taxonid=", species, "_model=", outacro, "_method=", model_name,".tif"),
                      overwrite = T)
          
          if (do_mess) {
            writeRaster(pred_mess, outmess, overwrite = T)
            cogeo_optim(outmess)
          }
          
          return(pred_f[["current"]])
          
        })
        
        treg <- obissdm::.get_time(treg, "Model prediction")
        
        
        
        # PART X: GET RESPONSE CURVES ----
        if (verbose) cli::cli_alert_info("Getting response curves")
        # Predict models
        model_respcurves <- lapply(good_models, function(id) {
          
          model_name <- model_fits[[id]]$name
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
          if (!dir.exists(outpath)) fs::dir_create(outpath)
          outfile <- paste0(outpath,
                            "taxonid=", species, "_model=", outacro, "_method=", model_name,
                            "_what=respcurves.parquet")
          
          rcurves <- resp_curves(model_fits[[id]],
                                 subset(env$layers, multi_mod_max$best_variables))
          
          write_parquet(rcurves, outfile)
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/figures/")
          if (!dir.exists(outpath)) fs::dir_create(outpath)
          outfile <- paste0(outpath,
                            "taxonid=", species, "_model=", outacro, "_method=", model_name,
                            "_what=responsecurves.png")
          
          p <- plot(rcurves) + ggtitle(paste("AphiaID:", species, "| Model:", outacro, "| Algo:", model_name))
          ggplot2::ggsave(filename = outfile, plot = p,
                          width = 10, height = 8)
          
          return(invisible(NULL))
          
        })
        
        treg <- obissdm::.get_time(treg, "Response curves")
        
        
        
        # PART X: ENSEMBLE ----
        if (length(good_models) > 1) {
          if (verbose) cli::cli_alert_info("Enough number of models, ensembling")
          
          to_ensemble <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/"),
                                    pattern = "\\.tif$", full.names = T)
          to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]
          to_ensemble <- lapply(to_ensemble, rast)
          lays <- nlyr(to_ensemble[[1]])
          
          ensemble <- lapply(1:lays, function(l){
            pull_l <- lapply(to_ensemble, function(r){
              r[[l]]
            })
            pull_l <- rast(pull_l)
            ensemble_models("median", pull_l)
          })
          
          ensemble_cv <- rast(lapply(ensemble, function(x) x[[2]]))
          ensemble <- rast(lapply(ensemble, function(x) x[[1]]))
          
          names(ensemble) <- names(ensemble_cv) <- names(to_ensemble[[1]])
          names(ensemble_cv) <- paste0(names(ensemble_cv), "_sd")
          
          outens <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/", 
                           "taxonid=", species, "_model=", outacro, "_method=ensemble.tif")
          
          writeRaster(c(ensemble, ensemble_cv), outens,
                      #wopt=list(verbose=TRUE),
                      overwrite = T)
          
          model_predictions <- c(model_predictions, ensemble[[1]])
          
          res_ens <- cogeo_optim(outens)
          
          # Evaluate ensemble
          ensemble_curr <- ensemble[[1]]
          
          fit_ens <- extract(ensemble_curr, sp_data$coord_training, ID = F)
          
          fit_ens_eval <- obissdm::eval_metrics(sp_data$training$presence, fit_ens[,1])
          
          if (!is.null(sp_data$eval_data)) {
            eval_ens <- extract(ensemble_curr, sp_data$coord_eval, ID = F)
            
            eval_ens_eval <- obissdm::eval_metrics(sp_data$eval_data$presence, eval_ens[,1])
            
            ensemble_metrics <- rbind(cbind(as.data.frame(t(fit_ens_eval)), origin = "fit_data"),
                                      cbind(as.data.frame(t(eval_ens_eval)), origin = "eval_data"))
          } else {
            ensemble_metrics <- cbind(as.data.frame(t(fit_ens_eval)), origin = "fit_data")
          }
          
          # Get average of models
          avg_fit_metrics <- as.data.frame(t(apply(do.call("rbind", lapply(good_models, function(id){
            t(apply(model_fits[[id]]$cv_metrics, 2, mean, na.rm = T))
          })), 2, mean, na.rm = T)))
          avg_fit_metrics$origin <- "avg_fit"
          
          if (!is.null(sp_data$eval_data)) {
            avg_eval_metrics <- as.data.frame(t(apply(do.call("rbind", lapply(good_models, function(id){
              t(model_fits[[id]]$eval_metrics)
            })), 2, mean, na.rm = T)))
            avg_eval_metrics$origin <- "avg_eval"
            avg_fit_metrics <- rbind(avg_fit_metrics, avg_eval_metrics)
          }
          
          ensemble_metrics <- rbind(ensemble_metrics, avg_fit_metrics)
          
          # TODO: add SD?
          
          write_parquet(ensemble_metrics,
                        paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/", 
                               "taxonid=", species, "_model=", outacro,
                               "_method=ensemble_what=metrics.parquet"))
        }
        
        # Ensemble response curves
        if (length(good_models) > 1) {
          to_ensemble <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/"),
                                    pattern = "what=resp", full.names = T)
          to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]
          to_ensemble <- lapply(to_ensemble, read_parquet)
          
          init_data <- to_ensemble[[1]]
          
          for (k in 2:length(to_ensemble)) {
            tdf <- data.frame(to_ensemble[[2]][,c("response")])
            colnames(tdf) <- paste0("response_", k)
            init_data <- cbind(init_data, tdf)
          }
          
          ens_resp <- apply(init_data[,startsWith(colnames(init_data), "resp")], 1, "median")
          
          ens_resp_curves <- cbind(init_data[,c("variable", "base")], response = ens_resp)
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
          outfile <- paste0(outpath,
                            "taxonid=", species, "_model=", outacro, "_method=ensemble",
                            "_what=respcurves.parquet")
          
          write_parquet(ens_resp_curves, outfile)
          
        }
        
        
        
        # PART X: SAVE BINARIZATION INFO ----
        if (verbose) cli::cli_alert_info("Saving binarization info")
        # Get thresholds
        # P10 threshold
        thresh_p10_mtp <- lapply(model_predictions, function(pred) {
          predv <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1,], ID = F)[,1]
          p10 <- ceiling(length(predv) * 0.9)
          data.frame(p10 = rev(sort(predv))[p10], mtp = min(na.omit(predv)))
        })
        
        # Max Specificity + Sensitivity
        thresh_maxss <- lapply(model_predictions, function(pred) {
          pred_p <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1,], ID = F)[,1]
          pred_a <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 0,], ID = F)[,1]
          
          pred_eval <- predicts::pa_evaluate(p = pred_p, a = pred_a)
          
          pred_thresh <- predicts::threshold(pred_eval)
          
          pred_thresh
        })
        
        if (length(model_predictions) > length(good_models)) {
          names(thresh_p10_mtp) <- names(thresh_maxss) <- c(algorithms[good_models], "ensemble")
        } else {
          names(thresh_p10_mtp) <- names(thresh_maxss) <- algorithms[good_models]
        }
        
        thresh_p10_mtp <- bind_rows(thresh_p10_mtp, .id = "model")
        thresh_maxss <- bind_rows(thresh_maxss, .id = "model")
        
        thresh <- left_join(thresh_p10_mtp, thresh_maxss, by = "model")
        
        outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
        write_parquet(thresh, paste0(outpath,
                                     "taxonid=", species, "_model=", outacro,
                                     "_what=thresholds.parquet"))
        
        
        
        # PART X: OPTIMIZE GEOTIFF (COG FORMAT) ----
        if (verbose) cli::cli_alert_info("Optimizing rasters (COG)")
        
        to_optimize <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions"), full.names = T)
        to_optimize <- to_optimize[!grepl("aux", to_optimize)]
        to_optimize <- to_optimize[!grepl("cog", to_optimize)]
        
        optimize_result <- lapply(to_optimize, cogeo_optim)
        optimize_result <- do.call("rbind", optimize_result)
        
        
        
        # PART X: CREATE MASKS AND SAVE ----
        if (verbose) cli::cli_alert_info("Saving masks")
        
        # Get base raster
        base_layer <- rast("data/env/terrain/bathymetry_mean.tif")
        base_layer[!is.na(base_layer)] <- 1
        
        # Get mask for study area species occurrence
        ecoreg_occ_mask <- mask(base_layer, ecoreg_occ)
        
        # Get mask for study area ecoregions
        ecoreg_mask <- mask(base_layer, ecoreg_sel)
        
        # Get area used for fitting
        fit_mask <- terra::extend(env$layers[[1]], model_predictions[[1]])
        fit_mask[!is.na(fit_mask)] <- 1
        
        masks <- c(ecoreg_occ_mask, ecoreg_mask, fit_mask)
        names(masks) <- c("native_ecoregions", "fit_ecoregions", "fit_region")
        
        outmask <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/", 
                          "taxonid=", species, "_model=", outacro, "_mask.tif")
        terra::writeRaster(masks, outmask, overwrite = T)
        mask_opt <- cogeo_optim(outmask)
        if (file.exists(paste0(outmask, ".aux.json"))) fs::file_delete(paste0(outmask, ".aux.json"))
        
        treg <- obissdm::.get_time(treg, "File optimization")
        
        
        
        # PART X: EVALUATE MODELS USING OTHER TECHNIQUES ----
        if (verbose) cli::cli_alert_info("Performing post-evaluation")
        
        # TEST 1: ENVIRONMENTAL SPACE
        # Make PCA
        # pca_env <- princomp(scale(subset(env$layers, multi_mod_max$best_variables)),
        #                     maxcell = (ncell(env$layers[[1]]) * .5))
        # pca_env_pred <- predict(scale(subset(env$layers, multi_mod_max$best_variables)), pca_env, index = 1:2)
        # 
        # # Extract at data points
        # pca_fit <- terra::extract(pca_env_pred, sp_data$coord_training[sp_data$training$presence == 1,], ID = F)
        # 
        # # Extract at SAMPLED data points
        # pca_pred <- lapply(model_predictions, function(pred){
        #   pred <- mask(pred, masks$fit_region)
        #   
        #   samp_pred <- spatSample(pred, size = sum(sp_data$training$presence), method = "weights", values = F, xy = T, na.rm = T)
        #   
        #   terra::extract(pca_env_pred, samp_pred, ID = F)
        # })
        # 
        # plot(NULL, xlim = minmax(pca_env_pred)[,1], ylim = minmax(pca_env_pred)[,2])
        # 
        # points(pca_fit, pch = 20, cex = .5, col = "blue")
        # points(pca_pred[[1]], pch = 20, cex = .5, col = "red")
        # 
        # chul1 <- chull(pca_fit)
        # chul1 <- c(chul1, chul1[1])
        # chul1 <- vect(as.matrix(pca_fit[chul1,]), type = "polygons")
        # 
        # chul2 <- chull(pca_pred[[1]])
        # chul2 <- c(chul2, chul2[1])
        # chul2 <- vect(as.matrix(pca_pred[[1]][chul2,]), type = "polygons")
        # 
        # lines(chul1, col = "blue")
        # lines(chul2, col = "red")
        
        
        # Add Hypervolume option, # TODO
        
        # TEST 2 - temperature kernel
        # Load temperature information
        temp_layer <- list.files("data/env/", recursive = T, 
                                 pattern = "thetao_baseline",
                                 full.names = T)
        temp_layer <- temp_layer[grepl(hab_depth, temp_layer)]
        temp_layer <- temp_layer[grepl("mean\\.tif$", temp_layer)]
        temp_layer <- rast(temp_layer)
        
        if (ext(temp_layer) != ext(model_predictions[[1]])) {
          temp_layer <- crop(temp_layer, model_predictions[[1]])
        }
        
        data_fit <- extract(temp_layer, 
                            dplyr::bind_rows(sp_data$coord_training,
                                             sp_data$coord_eval), ID = F)
        
        kd <- ks::kde(data_fit[,1], h = 8)
        percentiles <- ks::qkde(c(0.01, 0.99), kd)
        
        masked <- temp_layer
        masked[temp_layer >= percentiles[1] & temp_layer <= percentiles[2]] <- 3
        masked[temp_layer < percentiles[1] | temp_layer > percentiles[2]] <- 0
        
        binary_maps <- lapply(1:length(model_predictions), function(x){
          pred <- model_predictions[[x]]
          pred[pred < thresh_p10_mtp[x,]$p10] <- 0
          pred[pred > 0] <- 1
          pred
        })
        
        plot(masked)
        plot(binary_maps[[1]])
        
        diff_maps <- lapply(binary_maps, function(x){
          x + masked
        })
        
        diff_perc <- lapply(diff_maps, function(x){
          x[x == 0] <- NA
          x[x == 3] <- NA
          result <- as.data.frame(x)
          result <- as.data.frame(table(result[,1]))
          colnames(result) <- c("status", "frequency")
          result$status <- ifelse(result$status == 1, "outside_tenv", "inside_tenv")
          result$percentage <- round((result$frequency * 100) / sum(result$frequency), 2)
          result
        })
        
        trange_maps <- lapply(binary_maps, function(x){
          x[x != 1] <- NA
          x_pts <- as.data.frame(x, xy = T)[,1:2]
          quantile(extract(temp_layer, x_pts, ID = F)[,1],
                   c(0.01, 0.05, 0.5, 0.95, 0.99), na.rm = T)
        })
        
        for (gm in 1:length(good_models)) {
          salg <- algorithms[good_models[gm]]
          model_log$model_posteval[[salg]] <- list(
            thermal_range = quantile(data_fit[,1], c(0.01, 0.05, 0.5, 0.95, 0.99)),
            thermal_range_binary = trange_maps[[gm]],
            thermal_envelope = diff_perc[[gm]]
          )
        }
        
        # jacc <- function(raster1.bin, raster2.bin) {
        #   combination <- sum(raster1.bin, raster2.bin, na.rm= T)
        #   intersection <- combination == 2
        #   
        #   # Union is all the area covered by the both rasters
        #   union <- combination >= 1
        #   
        #   union <- as.data.frame(union)
        #   inter <- as.data.frame(intersection)
        #   
        #   length(inter[inter[,1],]) / length(union[union[,1],])
        # }
        # 
        # overlap_metric <- jacc(masked, crop(binary_maps[[1]], masked))
        # overlap_metric_b <- dismo::nicheOverlap(raster::raster(masked), raster::raster(crop(binary_maps[[1]], masked)))
        # 
        # 
        # Test 4: niche overlap ecospat
        
        
        
        # PART X: SAVE MODELS OBJECTS ----
        if (verbose) cli::cli_alert_info("Saving models objects")
        
        # Save metrics
        save_metrics <- lapply(good_models, function(id){
          model_name <- model_fits[[id]]$name
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
          outfile <- paste0(outpath,
                            "taxonid=", species, "_model=", outacro, "_method=", model_name)
          
          write_parquet(model_fits[[id]]$cv_metrics, paste0(outfile, "_what=cvmetrics.parquet"))
          
          other_metrics <-  list(model_fits[[id]]$full_metrics,
                                 NULL)#model_fits[[id]]$eval_metrics)
          names(other_metrics) <- c("full", "eval")
          other_metrics <- dplyr::bind_rows(other_metrics, .id = "what")
          
          write_parquet(other_metrics, paste0(outfile, "_what=fullmetrics.parquet"))
        })
        
        # Update log with best parameters
        for (al in 1:length(algorithms)) {
          if (!is.null(model_fits[[al]])) {
            model_log$model_bestparams[[al]] <- model_fits[[al]]$parameters
          }
        }
        
        # Save models
        save_models <- lapply(good_models, function(id) {
          
          model_name <- model_fits[[id]]$name
          
          outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
          if (!dir.exists(outpath)) fs::dir_create(outpath)
          outfile <- paste0(outpath,
                            "taxonid=", species, "_model=", outacro, "_method=", model_name,
                            "_what=model.rds")
          
          saveRDS(model_fits[[id]], file = outfile)
          
          return(invisible(NULL))
          
        })
        
        # Save fit points
        write_parquet(sp_data$coord_training[sp_data$training$presence == 1,],
                      paste0(outfolder, "/taxonid=", species, "/model=", outacro,
                             "/taxonid=", species, "_model=", outacro, "_",
                             "what=fitocc.parquet"))
        
        
        # Save log and return object
        outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/")
        outfile <- paste0(outpath,
                          "taxonid=", species, "_model=", outacro, "_what=log.json")
        save_log(model_log, outfile)
        
        cli::cli_alert_success("Model concluded for species {species}.")
        
        st_status <- "succeeded"
      } else {
        if (verbose) cli::cli_alert_warning("No good model available")
        st_status <- "no_good_model"
      }
      
    } else {
      if (verbose) cli::cli_alert_warning("All models failed")
      st_status <- "all_failed"
    }
    
  } else {
    if (nrow(species_data) < 15) {
      if (verbose) cli::cli_alert_warning("Low number of points, failed")
      # If data is low, annotate storr
      st_status <- "low_data"
    } else {
      if (verbose) cli::cli_alert_warning("No data available, failed")
      # If data is not available, annotate storr
      st_status <- "no_data"
    }
  }
  
  # Return nothing, results are saved
  return(st_status)
  
}
