############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################## Components for main modelling function ######################

# Loading main data ----
.cm_load_species <- function(species, species_dataset) {
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
  
  return(species_data)
}


.cm_check_coastal <- function(species_data, env) {
  
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
  
  return(env)
}


.cm_prepare_data_obj <- function(fit_pts, env, verbose = FALSE) {
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
  
  return(sp_data)
}



.cm_bias_assess <- function(sp_data, env) {
  
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
  
  return(list())
  
}



# .cm_fit_models <- function(sp_data, algorithms) {
#   
#   model_fits <- list()
#   for (i in 1:length(algorithms)) {
#     if (verbose) cli::cli_alert_info("Fitting {algorithms[i]}")
#     model_fits[[i]] <- try(do.call(paste0("sdm_module_", algorithms[i]), list(sp_data, verbose = verbose)))
#     if (inherits(model_fits[[i]], "try-error")) {
#       model_log$model_result[[algorithms[i]]] <- "failed"
#       model_fits[[i]] <- NULL
#     } else {
#       model_log$model_result[[algorithms[i]]] <- "succeeded"
#     }
#   }
#   
#   return(model_fits)
# }



.cm_predict_models <- function() {
  
  # Predict models
  model_predictions <- lapply(good_models, function(id) {
    model_name <- model_fits[[id]]$name
    
    best_hyp <- multi_mod_max$best_model
    
    outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/")
    
    if (!dir.exists(outpath)) fs::dir_create(outpath)
    
    outmess <- paste0(outpath, "taxonid=", species, "_model=", outacro, "_what=mess.tif")
    
    if (!file.exists(gsub("mess", "mess_cog", outmess))) {do_mess <- TRUE} else {do_mess <- FALSE} 
    
    scenarios <- data.frame(
      scenario = c("current", rep(c("ssp1", "ssp2", "ssp3", "ssp4", "ssp5"),
                                  each = 2)),
      year = c(NA, rep(c(2050, 2100), 5))
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
          as.data.frame(env_to_pred, xy = T)[to_mess,], 
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
  
  return(model_predictions)
}

.cm_get_respcurves <- function() {
  
}