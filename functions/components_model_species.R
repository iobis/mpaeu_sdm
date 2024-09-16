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
    lf <- lf[grepl(paste0("key=", species, ".parquet"), lf)]
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

# Check coastal ----
.cm_check_coastal <- function(species_data, env, coord_names, verbose) {
  
  if ("coastal" %in% names(env$hypothesis)) {
    # Test if is within europe/coastal
    europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
    europe_starea <- terra::ext(europe_starea)
    
    is_within <- terra::is.related(terra::vect(species_data[,coord_names], geom = coord_names),
                            europe_starea, "intersects")
    if (!all(is_within)) {
      env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
      env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
      if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
    } else {
      
      wave <- terra::rast("data/env/terrain/wavefetch.tif")
      wave_pts <- terra::extract(wave, species_data[,coord_names])
      
      if (sum(is.na(wave_pts[,2])) > ceiling(nrow(species_data) * .05)) {
        env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
        env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
        if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
      } else if ((nrow(species_data) - sum(is.na(wave_pts[,2]))) < 15) {
        env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
        env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
        if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
      } else {
        env$layers <- terra::crop(env$layers, europe_starea)
      }
    }
  }
  
  return(env)
}

# Check ecoregions ----
.cm_check_ecoregions <- function(ecoregions, fit_pts, eval_pts,
                                 env, limit_by_depth, depth_buffer) {
  
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
    
  } else {
    if ("coastal" %in% names(env$hypothesis)) {
      ecoreg_sel <- crop(ecoreg_sel, europe_starea)
      env$layers <- mask(env$layers, ecoreg_sel)
    } else {
      env$layers <- mask(env$layers, ecoreg_sel)
    }
  }
  
  return(env)
}

# Prepare data object ----
.cm_prepare_data_obj <- function(fit_pts, eval_pts, env, verbose = FALSE) {
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


# Assess bias (TODO) ----
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

# Check good models ----
.cm_check_good_models <- function(model_fits, tg_metrics = "cbi", tg_threshold = 0.3) {
  
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
  
  return(good_models)
}

# Predict models ----
.cm_predict_models <- function(good_models, model_fits, best_hyp, hab_depth,
                               outfolder, species, outacro,
                               verbose) {
  
  # Predict models
  model_predictions <- lapply(good_models, function(id) {
    model_name <- model_fits[[id]]$name
    
    #best_hyp <- multi_mod_max$best_model
    
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
  
  return(model_predictions)
}

.cm_get_respcurves <- function() {
  
}

# Post evaluation ----
.cm_posteval_sst <- function(model_predictions, sdm_data,
                             thresholds, algorithms, hab_depth,
                             good_models, model_log) {
  
  temp_layer <- list.files("data/env/", recursive = T, 
                           pattern = "thetao_baseline",
                           full.names = T)
  temp_layer <- temp_layer[grepl(hab_depth, temp_layer)]
  temp_layer <- temp_layer[grepl("mean\\.tif$", temp_layer)]
  temp_layer <- terra::rast(temp_layer)
  
  if (ext(temp_layer) != ext(model_predictions[[1]])) {
    temp_layer <- terra::crop(temp_layer, model_predictions[[1]])
    temp_layer <- terra::mask(temp_layer, model_predictions[[1]])
  }
  
  data_fit <- terra::extract(temp_layer, 
                      dplyr::bind_rows(sdm_data$coord_training,
                                       sdm_data$coord_eval), ID = F)
  
  kd <- ks::kde(data_fit[,1], h = 8)
  percentiles <- ks::qkde(c(0.01, 0.99), kd)
  
  masked <- temp_layer
  masked[temp_layer >= percentiles[1] & temp_layer <= percentiles[2]] <- 3
  masked[temp_layer < percentiles[1] | temp_layer > percentiles[2]] <- 0
  
  binary_maps <- lapply(1:length(model_predictions), function(x){
    pred <- model_predictions[[x]]
    pred[pred < thresholds[x,]$p10] <- 0
    pred[pred > 0] <- 1
    pred
  })
  
  # plot(masked)
  # plot(binary_maps[[1]])
  
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
  
  for (gm in 1:length(model_predictions)) {
    salg <- c(algorithms[good_models], "ensemble")[gm]
    model_log$model_posteval[[salg]] <- list(
      thermal_range = quantile(data_fit[,1], c(0.01, 0.05, 0.5, 0.95, 0.99)),
      thermal_range_binary = trange_maps[[gm]],
      thermal_envelope = diff_perc[[gm]]
    )
    names(model_log$model_posteval[[salg]]$thermal_range) <- c("q_0.01", "q_0.05", "q_0.5", "q_0.95", "q_0.99")
    names(model_log$model_posteval[[salg]]$thermal_range_binary) <- c("q_0.01", "q_0.05", "q_0.5", "q_0.95", "q_0.99")
  }
  
  return(model_log)
}

.cm_posteval_hypervolume <- function(idx, var_imp, model_predictions,
                                     env_layers, sdm_data, n_vars = 3,
                                     return_plot = TRUE) {
  
  sink(nullfile())
  require(hypervolume)
  require(ggplot2)
  require(patchwork)
  
  # Select N most important vars
  vars <- var_imp[[idx]]
  
  if (n_vars > nrow(vars)) {
    n_vars <- nrow(vars)
  } else if (n_vars < 2) {
    n_vars < 2
  }
  
  vars <- vars$variable[1:n_vars]
  
  # Subset environmental layers
  climate <- terra::subset(env_layers, vars)
  climate <- terra::scale(climate)
  
  # Create hypervolumes
  data_o <- terra::extract(climate, sdm_data$coord_training[sdm_data$training$presence == 1,], ID = F)
  hyper_list <- list()
  hyper_list[[1]] <- suppressMessages(hypervolume(data_o, method = "gaussian", name = "original", verbose = FALSE))
  
  # Sample new occurrences
  # Remove very low values to avoid points on very unlikely areas
  predicted_dist <- model_predictions[[idx]]
  predicted_dist[predicted_dist < 0.05] <- 0
  
  pred_occ <- terra::spatSample(
    terra::mask(terra::crop(predicted_dist, env_layers[[1]]),
                env_layers[[1]]),
    size = sum(sdm_data$training$presence),
    method = "weights", na.rm = T, xy = T)
  
  data_p <- terra::extract(climate, pred_occ[,1:2], ID = F)
  
  hyper_list[[2]] <- suppressMessages(hypervolume(data_p, method = "gaussian", name = "predicted", verbose = FALSE))
  
  # Transform hyper_list in an HypervolumeList
  hyper_list <- hypervolume_join(hyper_list)
  
  # Get occupancy statistics
  hyper_occupancy <- hypervolume_n_occupancy(hyper_list, 
                                             classification = c("original", "predicted"),
                                             method = "box",
                                             FUN = "mean", verbose = F)
  
  climatic_occupancy <- hypervolume_to_data_frame(hyper_occupancy) %>%
    filter(ValueAtRandomPoints != 0) %>%
    group_by(Name) %>%
    sample_n(500) %>%
    ungroup() %>%
    sample_n(nrow(.)) %>%
    rename(type = Name, occupancy = ValueAtRandomPoints)
  
  sink()
  # Get plots
  if (return_plot) {
    
    plot_list <- list()
    
    var_comb <- expand.grid(option1 = colnames(climatic_occupancy)[2:(ncol(climatic_occupancy)-1)],
                            option2 = colnames(climatic_occupancy)[2:(ncol(climatic_occupancy)-1)])
    var_comb <- var_comb[var_comb$option1 != var_comb$option2, ]
    var_comb <- var_comb[!duplicated(t(apply(var_comb, 1, sort))), ]
    
    for (pt in 1:nrow(var_comb)) {
      
      eval(parse(text = paste0('
      c_hull <- na.omit(climatic_occupancy) %>%
        group_by(type) %>%
        mutate(hull = 1:n(), hull = factor(hull, chull(', var_comb[pt,1],', ', var_comb[pt,2],'))) %>%
        arrange(hull)
                               ')))
      
      eval(parse(text = paste0('plot_list[[pt]] <- ggplot(climatic_occupancy) +
        geom_point(aes(',var_comb[pt,1],', ',var_comb[pt,2],', color = type), alpha = 0.5) +
        geom_polygon(data = filter(c_hull, !is.na(hull)), aes(',var_comb[pt,1],', ',var_comb[pt,2],',
        fill = type), alpha = 0.1) +
        theme_bw() +
        scale_size_continuous(limits = c(0, 1)) +
        scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
        scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
        labs(title = "',paste0(var_comb[pt,1],' - ', var_comb[pt,2]),'")')))
    }
    
    p <- eval(parse(text = paste0(
      paste0("plot_list[[", 1:length(plot_list), "]]", collapse = "+"), "+ patchwork::plot_layout(guides = 'collect')"
    )))
  } else {
    p <- NULL
  }
  
  sink(nullfile())
  hv_set <- hypervolume_set(hyper_list[[1]], hyper_list[[2]], check.memory=FALSE, verbose = F)
  
  overlap <- hypervolume_overlap_statistics(hv_set)
  sink()
  
  return(list(
    overlap = overlap,
    occupancy = climatic_occupancy,
    plot = p
  ))
  
}

.cm_posteval_nicheequiv <- function(sdm_data, predicted_dist, env_layers,
                                    iterations = 5, plot_example = F,
                                    ...){
  #tictoc::tic()
  bckg_points <- sdm_data$training[sdm_data$training$presence == 0,-1]
  # Reduce if very large background sample
  if (nrow(bckg_points) > 50000) {
    bckg_points <- bckg_points[sample(1:nrow(bckg_points), 50000),]
  } else if (nrow(bckg_points) < 10000) {
    warning("Low number of background points for niche ecospat. Sample retrieved.")
    env_s <- terra::subset(env_layers, colnames(sdm_data$training)[-1])
    env_s <- as.data.frame(env_s)
    env_s <- na.omit(env_s)
    bckg_points <- env_s[sample(1:nrow(env_s), 10000),]
  }
  
  pca_env <- ade4::dudi.pca(bckg_points, 
                           center = T, scale = T, scannf = F, nf = 2)
  
  # Extract scores for background, in this case = env space
  scores_bckg <- pca_env$li
  
  # Extract scores occurrences
  scores_occ <- ade4::suprow(pca_env, sdm_data$training[sdm_data$training$presence == 1,-1])$lisup
  
  results <- data.frame(matrix(ncol = 10, nrow = iterations))
  colnames(results) <- c("D", "I", "expansion", "stability", "unfiling",
                         "p_D", "p_I", "p_expansion", "p_stability", "p_unfiling")
  
  # Remove very low values to avoid points on very unlikely areas
  predicted_dist[predicted_dist < 0.05] <- 0
  
  for (i in 1:iterations) {
    pred_occ <- terra::spatSample(
      terra::mask(terra::crop(predicted_dist, env_layers[[1]]), env_layers[[1]]),
      size = sum(sdm_data$training$presence),
      method = "weights", na.rm = T, xy = T)
    
    # Extract scores predicted
    scores_pred <- ade4::suprow(pca_env, terra::extract(terra::subset(env_layers, 
                                                                      colnames(sdm_data$training[,-1])),
                                                        pred_occ[,1:2], ID = F))$lisup
    
    # Calculate occurence density
    real_sp <- ecospat::ecospat.grid.clim.dyn(scores_bckg, scores_bckg, scores_occ, R = 100, ...)
    pred_sp <- ecospat::ecospat.grid.clim.dyn(scores_bckg, scores_bckg, scores_pred, R = 100, ...)
    
    # Test equivalency
    equivalency <- ecospat::ecospat.niche.equivalency.test(real_sp, pred_sp, rep = 10)
    
    results[i,] <- c(unlist(equivalency$obs), equivalency$p.D, equivalency$p.I,
                     equivalency$p.expansion, equivalency$p.stability, equivalency$p.unfilling)
  }
  #tictoc::toc()
  if (plot_example) ecospat::ecospat.plot.niche.dyn(real_sp, pred_sp)
  
  return(round(results, 3))
}


# Auxiliary functions ----
# Convert eco info to habitat depth
hab_to_depth <- function(hab, default = "depthsurf") {
  
  if (grepl("benthic|demersal|bottom", hab)) {
    dstrata <- "depthmean"
  } else if (grepl("pelagic_surface", hab)) {
    dstrata <- "depthsurf"
  } else if (grepl("NOT_FOUND", hab)) {
    dstrata <- default
  } else {
    dstrata <- "depthsurf"
  }
  
  return(dstrata)
}

# Split dataset into eval and fit
split_ds <- function(sp_occ,
                     what = "fit",
                     only_coords = TRUE) {
  
  if (grepl("fit", what)) {
    to_get <- "fit_points"
  } else if (grepl("eval", what)) {
    to_get <- "eval_points"
  } else {
    stop("What should be one of `fit_points` or `eval_points`")
  }
  
  pts <- sp_occ %>%
    filter(data_type == to_get) %>%
    as.data.frame()
  
  if (nrow(pts) == 0) {
    pts <- NULL
  } else if (only_coords) {
    pts <- pts[,c("decimalLongitude", "decimalLatitude")]
  }
  
  return(pts)
}

# Reduce file for saving
.reduce_model_file <- function(sdm) {
  if (grepl("gam_", sdm$name)) {
    sdm$model$model <- NULL
    sdm$model$wt <- NULL
    sdm$model$y <- NULL
    sdm$model$prior.weights <- NULL
    sdm$model$residuals <- NULL
    sdm$model$fitted.values <- NULL
    sdm$model$linear.predictors <- NULL
    sdm$model$weights <- NULL
    sdm$model$offset <- NULL
  } else if (grepl("rf_", sdm$name)) {
    #From https://stats.stackexchange.com/questions/102667/reduce-random-forest-model-memory-size
    # This part applies only to caret models
    sdm$model$finalModel$predicted <- NULL 
    sdm$model$finalModel$oob.times <- NULL 
    sdm$model$finalModel$y <- NULL
    sdm$model$finalModel$votes <- NULL
    sdm$model$control$indexOut <- NULL
    sdm$model$control$index    <- NULL
    sdm$model$trainingData <- NULL
    # This applies to all
    sdm$model$err.rate <- NULL
    sdm$model$predicted <- NULL
    sdm$model$votes <- NULL
    sdm$model$oob.times <- NULL
    sdm$model$y <- NULL
    attr(sdm$model$terms,".Environment") <- c()
    attr(sdm$model$formula,".Environment") <- c()
  } else if (sdm$name == "maxent") {
    sdm$model$dev.ratio <- NULL
    sdm$model$beta <- NULL
  } else if (sdm$name == "xgboost") {
    sdm$model$evaluation_log <- NULL
    sdm$model$raw <- NULL
  }
  
  return(sdm)
}

# Normalize raster
# .normalize_raster <- function(pred, range = NULL) {
#   if (is.null(range)) {
#     range <- terra::minmax(pred)[,1]
#   }
#   pred <- (pred - min(range)) / (max(range) - min(range))
#   return(list(prediction = pred,
#               range = range))
# }
