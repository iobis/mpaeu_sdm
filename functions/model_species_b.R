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
                          limit_by_depth = TRUE,
                          depth_buffer = 500,
                          assess_bias = TRUE,
                          correct_bias = TRUE,
                          verbose = FALSE) {
  
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
  
  # PART 1: DATA LOADING ----
  # Record timing
  treg <- obissdm::.get_time()
  
  # Define global objects and fill log
  coord_names <- c("decimalLongitude", "decimalLatitude")
  model_log <- obissdm::gen_log(algorithms)
  model_log$taxonID <- species
  model_log$model_acro <- outacro
  model_log$model_date <- Sys.Date()
  
  # Load species data
  if (verb_1) cli::cli_alert_info("Reading data")
  species_data <- .cm_load_species(species, species_dataset)
  
  treg <- obissdm::.get_time(treg, "Species data loading")
  
  
  
  # If data is available, proceeds
  if (nrow(species_data) >= 15) {
    
    if (verb_1) cli::cli_alert_info("Enough number of points, proceeding")
    # Load species information
    species_name <- species_data$species[1]
    model_log$scientificName <- species_name
    
    # Load ecological information
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
    
    # Load environmental data
    env <- get_envofgroup(group, depth = hab_depth, load_all = T,
                          conf_file = "sdm_conf.yml", verbose = verbose)
    
    env <- .cm_check_coastal(species_data, env, coord_names, verb_1)
    
    # Load ecoregions shapefile
    ecoregions <- vect("data/shapefiles/MarineRealms.shp")
    
    # Split the dataset
    fit_pts <- split_ds(species_data)
    eval_pts <- split_ds(species_data, "eval")
    
    treg <- obissdm::.get_time(treg, "Data loading")
    
    
    
    # PART 2: DATA PREPARING ----
    if (verbose) cli::cli_alert_info("Preparing data")
    
    # Ecoregions checking
    env <- .cm_check_ecoregions(
      ecoregions, fit_pts, eval_pts, env,
      limit_by_depth, depth_buffer
    )
    
    model_log$model_details$limited_by_depth <- limit_by_depth
    model_log$model_details$depth_buffer <- depth_buffer
    
    # Prepare data object
    sp_data <- .cm_prepare_data_obj(fit_pts, eval_pts = eval_pts,
                                    env, verbose = verb_2)
    
    model_log$model_details$block_size <- unname(sp_data$blocks$grid_resolution)
    model_log$model_fit_points <- sum(sp_data$training$presence)
    model_log$model_eval_points <- sum(sp_data$eval_data$presence)
    
    treg <- obissdm::.get_time(treg, "Data preparing")
    
    
    
    # PART 3: ASSESS SPATIAL BIAS ----
    if (assess_bias) {
      if (verb_1) cli::cli_alert_info("Assessing spatial bias")
      
      # Add some way of assessing spatial bias
      model_log$model_details$control_bias <- correct_bias
      
      bias_todo <- .cm_bias_assess(sp_data, env)
        
      treg <- obissdm::.get_time(treg, "Sample bias assessment")
    }
    
    
    
    # PART 4: MODEL SELECTION ----
    if (verb_1) cli::cli_alert_info("Model selection")
    multi_mod_max <- sdm_multhypo(sdm_data = sp_data, sdm_method = "maxent",
                                  variables = env$hypothesis,
                                  verbose = verb_2)
    
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
    if (verb_1) cli::cli_alert_info("Starting model fitting")
    model_fits <- list()
    for (i in 1:length(algorithms)) {
      if (verbose) cli::cli_alert_info("Fitting {algorithms[i]}")
      model_fits[[i]] <- try(do.call(paste0("sdm_module_", algorithms[i]), list(sp_data, verbose = verb_2)))
      if (inherits(model_fits[[i]], "try-error")) {
        model_log$model_result[[algorithms[i]]] <- "failed"
        model_fits[[i]] <- NULL
      } else {
        model_log$model_result[[algorithms[i]]] <- "succeeded"
      }
    }
    
    treg <- obissdm::.get_time(treg, "Model fit")
    
    
    
    if (any(unlist(model_log$model_result) == "succeeded")) {
      if (verb_1) cli::cli_alert_info("Model fitting concluded. Checking if models are good")
      # PART 5: PREDICTION ----
      
      good_models <- .cm_check_good_models(model_fits, 
                                           tg_metrics = "cbi",
                                           tg_threshold = 0.3)
      
      if (length(good_models) > 0) {
        if (verbose) cli::cli_alert_info("Good models available, predicting.")
        
        # Predict models
        model_predictions <- .cm_predict_models(
          good_models, model_fits,
          best_hyp = multi_mod_max$best_model, hab_depth,
          outfolder, species, outacro,
          verbose = verb_2
        )
        
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
