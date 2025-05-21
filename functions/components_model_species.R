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
    # to do!
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

# Check habitat
.cm_check_habitat <- function(species) {
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

  return(eco_info$mode_life)
}



# Check coastal ----
.cm_check_coastal <- function(species_data, env, coord_names, verbose) {
  
  if ("coastal" %in% names(env$hypothesis)) {
    # # Test if is within europe/coastal
    # europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
    # europe_starea <- terra::ext(europe_starea)
    
    # is_within <- terra::is.related(terra::vect(species_data[,coord_names], geom = coord_names),
    #                         europe_starea, "intersects")
    # if (!all(is_within)) {
    #   env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
    #   env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
    #   if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
    # } else {
      
    #   wave <- terra::rast("data/env/terrain/wavefetch.tif")
    #   wave_pts <- terra::extract(wave, species_data[,coord_names])
      
    #   if (sum(is.na(wave_pts[,2])) > ceiling(nrow(species_data) * .05)) {
    #     env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
    #     env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
    #     if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
    #   } else if ((nrow(species_data) - sum(is.na(wave_pts[,2]))) < 15) {
    #     env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
    #     env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
    #     if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded")
    #   } else {
    #     env$layers <- terra::crop(env$layers, europe_starea)
    #   }
    # }
    wave <- terra::rast("data/env/terrain/wavefetch.tif")
    wave_pts <- terra::extract(wave, species_data[,coord_names])
    
    if (sum(is.na(wave_pts[,2])) > ceiling(nrow(species_data) * .05)) {
      env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
      env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
      if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded. Too many points out of area.")
    } else if ((nrow(species_data) - sum(is.na(wave_pts[,2]))) < 15) {
      env$hypothesis <- env$hypothesis[names(env$hypothesis) != "coastal"]
      env$layers <- terra::subset(env$layers, "wavefetch", negate = T)
      if (verbose) cli::cli_alert_warning("Coastal hypothesis found but discarded. Final points would be less than 15.")
    }
  }
  
  return(env)
}

# Check ecoregions ----
.cm_check_ecoregions <- function(ecoregions, fit_pts, eval_pts,
                                 env, limit_by_depth, depth_buffer,
                                 coord_names, model_log, verb_1) {
  
  # Check which eco-regions are covered by the points
  ecoreg_unique <- ecoregions$Realm[is.related(ecoregions,
                                               terra::vect(bind_rows(fit_pts, eval_pts), 
                                                    geom = coord_names), 
                                               "intersects")]
  
  model_log$model_details$ecoregions <- ecoreg_unique
  ecoreg_occ <- ecoregions[ecoregions$Realm %in% ecoreg_unique,]
  
  # Apply a buffer to ensure that all areas are covered
  sf::sf_use_s2(FALSE)
  ecoreg_occ_buff <- suppressMessages(
    suppressWarnings(terra::vect(sf::st_buffer(sf::st_as_sf(vect(bind_rows(fit_pts, eval_pts),
                                                          geom = coord_names)), 0.2)))
  )

  adj_ecoreg <- ecoregions$Realm[is.related(ecoregions, ecoreg_occ_buff,
                                            "intersects")]
  
  # Mask areas
  ecoreg_sel <- ecoregions[ecoregions$Realm %in% unique(
    c(ecoreg_unique, adj_ecoreg)
  ),]
  model_log$model_details$ecoregions_included <- unique(ecoreg_sel$Realm)
  
  ecoreg_sel <- terra::crop(ecoreg_sel, terra::ext(-180, 180, -90, 90))

  # Load bathymetry layer
  bath <- terra::rast("data/env/terrain/bathymetry_mean.tif")
  bath <- terra::mask(bath, ecoreg_sel)
  
  bath_pts <- terra::extract(bath, bind_rows(fit_pts, eval_pts))
  bath_pts <- bath_pts[!is.na(bath_pts[,2]),]
  
  # Limit by depth if TRUE
  if (limit_by_depth) {
    if (verb_1) cli::cli_alert_info("Limiting by depth")
    
    bath_range <- range(bath_pts[,2])
    bath_range[1] <- bath_range[1] - depth_buffer
    bath_range[2] <- ifelse((bath_range[2] + depth_buffer) > 0, 
                            0, bath_range[2] + depth_buffer)
    
    bath[bath < bath_range[1] | bath > bath_range[2]] <- NA
    
    if ("coastal" %in% names(env$hypothesis)) {
      #europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
      #bath <- terra::crop(bath, europe_starea)
      env$layers <- terra::mask(env$layers, bath)
      env$layers <- terra::mask(env$layers, env$layers$wavefetch)
    } else {
      env$layers <- terra::mask(env$layers, bath)
    }
    
    model_log$model_details$limited_by_depth <- TRUE
    model_log$model_details$depth_buffer <- depth_buffer
    
  } else {
    if ("coastal" %in% names(env$hypothesis)) {
     # europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
      #ecoreg_sel <- terra::crop(ecoreg_sel, europe_starea)
      env$layers <- terra::mask(env$layers, ecoreg_sel)
      env$layers <- terra::mask(env$layers, env$layers$wavefetch)
    } else {
      env$layers <- terra::mask(env$layers, ecoreg_sel)
    }
  }

  model_log$range_depth <- range(bath_pts[,2])
  
  return(list(
    env = env,
    model_log = model_log,
    bath_pts = bath_pts,
    ecoreg_occ = ecoreg_occ,
    ecoreg_sel = ecoreg_sel
  ))
}

# Check shallow ----
.cm_check_shallow <- function(bath_pts, env, verb_1) {
  if (min(bath_pts[,2], na.rm = T) >= -50) {
    if (verb_1) cli::cli_alert_info("Species from shallow areas - removing bathymetry/distance to coast")
    if (any(grepl("bathymetry", names(env$layers)))) {
      env$layers <- terra::subset(env$layers, "bathymetry_mean", negate = T)
      env$hypothesis <- lapply(env$hypothesis, function(x) x[!grepl("bathymetry", x)])
    }
    if (any(grepl("distcoast", names(env$layers)))) {
      env$layers <- terra::subset(env$layers, "distcoast", negate = T)
      env$hypothesis <- lapply(env$hypothesis, function(x) x[!grepl("distcoast", x)])
    }
  }
  return(env)
}

# Calculate quad samp -----
.cm_calc_quad <- function(env, quad_samp, fit_pts_sac) {
  env_size_t <- nrow(as.data.frame(env$layers, xy = T, na.rm = T))
  if (quad_samp <= 1 && quad_samp > 0) {
    quad_samp <- ceiling(env_size_t * quad_samp)
    if (quad_samp < 10000) quad_samp <- 10000
  }
  # Check if size of dataset is smaller then quad_samp
  if (nrow(fit_pts_sac) >= quad_samp) {
    est_tsize <- nrow(fit_pts_sac) * 2
    if (est_tsize > env_size_t) {
      quad_samp <- env_size_t
      est_tsize <- floor(quad_samp * 0.9)
      fit_pts_sac <- fit_pts_sac %>% dplyr::slice_sample(n = est_size)
      model_log$other_details <- c(model_log$other_details,
                                    paste("fit points sampled to", est_tsize))
    } else {
      quad_samp <- est_tsize
    }
  }
  quad_n <- ifelse(env_size_t < quad_samp, round(env_size_t), quad_samp)
  return(quad_n)
}


# Assess spatial bias
.cm_spat_bias <- function(metric_out, species, outacro, env, sp_data) {
  require(spatstat)
      
  # Get spatstat statistics and point process for control
  spat_im <- as.im(as.data.frame(aggregate(env$layers[[1]], 10, na.rm = T), xy = T))
  spat_window <- as.owin(spat_im)
  
  # Define ppp object based on point locations
  sppcords <- sp_data$coord_training[sp_data$training$presence == 1,]
  
  if (nrow(sppcords) > 10000) {
    sppcords <- sppcords %>% dplyr::slice_sample(n = 10000)
  }
  
  spat_ppp <- ppp(x = sppcords$decimalLongitude,
                  y = sppcords$decimalLatitude,
                  window = spat_window)
  
  # calculate envelope around L-hat estimates.
  spat_env_k <- envelope(spat_ppp, Kest, verbose = F)
  spat_env_l <- envelope(spat_ppp, Lest, verbose = F)
  
  # save information
  if (!dir.exists(metric_out)) fs::dir_create(metric_out)
  outfile <- paste0(metric_out, "taxonid=", species, "_model=", outacro,
                    "_what=biasmetrics.rds")
  saveRDS(list(k_stat = spat_env_k, l_stat = spat_env_l),
          file = outfile)

  return(invisible(NULL))
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


# Fit model -----
.cm_model_fit <- function(algorithms, algo_opts, sp_data,
                          model_log, verb_1, verb_2) {
  if (!is.null(algo_opts)) {
    algo_opts <- algo_opts[algorithms]
  } else {
    algo_opts <- obissdm::sdm_options()
    algo_opts <- algo_opts[algorithms]
  }
  model_fits <- lapply(seq_along(algorithms), \(x) NULL)
  for (i in seq_len(length(algorithms))) {
    if (verb_1) cli::cli_alert_info("Fitting {algorithms[i]}")
    model_fits[[i]] <- try(do.call(paste0("sdm_module_", algorithms[i]),
                                    list(sp_data, verbose = verb_2,
                                        options = algo_opts[[i]])))
    if (inherits(model_fits[[i]], "try-error")) {
      model_log$model_result[[algorithms[i]]] <- "failed"
      model_fits[[i]] <- FALSE
    } else {
      model_log$model_result[[algorithms[i]]] <- "succeeded"
    }
  }

  return(list(
    fits = model_fits,
    logs = model_log
  ))
}



# Check good models ----
.cm_check_good_models <- function(model_fits, tg_metric, tg_threshold) {

  # Get what models are above threshold
  good_models <- lapply(model_fits, function(model){
    if (!is.logical(model)) {
      cv_res <- model$cv_metrics
      cv_res <- apply(cv_res, 2, mean, na.rm = T)
      the_metric <- cv_res[[tg_metric]]
      if (!is.na(the_metric)) {
        if (the_metric >= tg_threshold) {
          if (sum(is.na(model$cv_metrics[[tg_metric]])) >= 4) {
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
.cm_predict_models <- function(good_models, model_fits, multi_mod,
                               pred_out, group, hab_depth,
                               sp_data, outfolder, outacro, species, env,
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

    scenarios <- data.frame(
      scenario = c("current", rep(c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
        each = 2
      )),
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
      outpred <- paste0(
        pred_out, "taxonid=", species, "_model=", outacro,
        "_method=", model_name, "_scen=", scenarios$scenario[k],
        ifelse(is.null(period), "", paste0("_", period)), ".tif"
      )

      env_to_pred <- obissdm::get_envofgroup(group,
        depth = hab_depth, load_all = F,
        scenario = scenarios$scenario[k],
        period = period,
        hypothesis = best_hyp,
        env_folder = "data/env",
        conf_file = "sdm_conf.yml",
        verbose = verb_2
      )
      env_to_pred <- terra::subset(env_to_pred, colnames(sp_data$training)[-1])

      if (best_hyp == "coastal") {
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

.cm_get_respcurves <- function(good_models, model_fits, sp_data,
                               env, multi_mod, metric_out, fig_out,
                               species, outacro) {
  lapply(good_models, function(id) {
          
    model_name <- model_fits[[id]]$name
    
    if (!dir.exists(metric_out)) fs::dir_create(metric_out)
    outfile <- paste0(metric_out,
                      "taxonid=", species, "_model=", outacro, "_method=", model_name,
                      "_what=respcurves.parquet")
    
    rcurves <- resp_curves(model_fits[[id]],
                            subset(env$layers, multi_mod$best_variables),
                            sdm_data = sp_data)
    
    arrow::write_parquet(rcurves, outfile)
    
    if (!dir.exists(fig_out)) fs::dir_create(fig_out)
    outfile <- paste0(fig_out,
                      "taxonid=", species, "_model=", outacro, "_method=", model_name,
                      "_what=responsecurves.png")
    
    p <- plot(rcurves) + ggtitle(paste("AphiaID:", species, "| Model:", outacro, "| Algo:", model_name))
    ggplot2::ggsave(filename = outfile, plot = p,
                    width = 10, height = 8)
    
    return(invisible(NULL))
    
  })
}

.cm_get_importance <- function(good_models, model_fits, 
                               sp_data, metric_out,
                               species, outacro) {
  model_varimport <- lapply(good_models, function(id){
    model_name <- model_fits[[id]]$name
    varimport <- obissdm::variable_importance(model_fits[[id]], sdm_data = sp_data)
    
    outfile <- paste0(metric_out,
                      "taxonid=", species, "_model=", outacro, "_method=", model_name,
                      "_what=varimportance.parquet")
    
    arrow::write_parquet(varimport, outfile)
    return(varimport)
  })

  if (length(good_models) > 1) {
    model_varimport_means <- lapply(model_varimport, function(x) {
      x$mean
    })
    model_varimport_sd <- apply(do.call("cbind", model_varimport_means), 
                                1, sd, na.rm = T)
    model_varimport_sd <- round(model_varimport_sd, 3)
    model_varimport_means <- apply(do.call("cbind", model_varimport_means), 
                                    1, mean, na.rm = T)
    model_varimport_means <- round(model_varimport_means, 3)
    ens_var_imp <- model_varimport[[1]][,1:2]
    ens_var_imp$mean <- model_varimport_means
    ens_var_imp$sd <- model_varimport_sd
    outfile <- paste0(metric_out,
                      "taxonid=", species, "_model=", outacro, "_method=ensemble",
                      "_what=varimportance.parquet")
    
    arrow::write_parquet(ens_var_imp, outfile)
    model_varimport <- c(model_varimport, list(ens_var_imp))
  }

  return(model_varimport)
}

.cm_ensemble_models <- function(species, good_models, algorithms, outfolder, 
                                sp_data, model_fits, metric_out,
                                pred_out, outacro, model_log, model_predictions, verb_1) {
  if (verb_1) cli::cli_alert_info("Enough number of models, ensembling")

  to_ensemble <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/"),
    pattern = "\\.tif$", full.names = T
  )
  to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]

  to_ensemble_g <- stringr::str_extract(to_ensemble, "scen=*.*_")
  to_ensemble_g <- gsub("scen=", "", to_ensemble_g)
  to_ensemble_g <- unique(to_ensemble_g)

  for (un in seq_len(length(to_ensemble_g))) {
    to_ensemble_r <- rast(to_ensemble[grepl(to_ensemble_g[un], to_ensemble)])

    ensemble <- ensemble_models("median", to_ensemble_r)

    ensemble_cv <- ensemble[[2]]
    ensemble_mean <- ensemble[[1]]

    if (grepl("current", to_ensemble_g[un])) {
      ensemble_curr <- ensemble_mean
    }

    names(ensemble_mean) <- gsub("_$", "", to_ensemble_g[un])
    names(ensemble_cv) <- paste0(names(ensemble_mean), "_sd")

    outens <- paste0(
      pred_out,
      "taxonid=", species, "_model=", outacro, "_method=ensemble_scen=", gsub("_$", "", to_ensemble_g[un]), ".tif"
    )

    ensemble_mean <- terra::as.int(ensemble_mean)
    ensemble_cv <- terra::as.int(ensemble_cv)
    terra::writeRaster(c(ensemble_mean, ensemble_cv), outens,
      overwrite = T, datatype = "INT1U"
    )
    res_ens <- cogeo_optim(outens)
    rm(ensemble, ensemble_mean, ensemble_cv)
  }

  model_predictions <- c(model_predictions, ensemble_curr)
  model_log$model_result$ensemble <- "succeeded"

  # Evaluate ensemble
  ensemble_curr <- ensemble_curr / 100

  fit_ens <- terra::extract(ensemble_curr, sp_data$coord_training, ID = F)

  fit_ens_eval <- obissdm::eval_metrics(sp_data$training$presence, fit_ens[, 1])

  if (!is.null(sp_data$eval_data)) {
    eval_ens <- terra::extract(ensemble_curr, sp_data$coord_eval, ID = F)

    eval_ens_eval <- obissdm::eval_metrics(sp_data$eval_data$presence, eval_ens[, 1])

    ensemble_metrics <- rbind(
      cbind(as.data.frame(t(fit_ens_eval)), what = "full", origin = "fit_data"),
      cbind(as.data.frame(t(eval_ens_eval)), what = "full", origin = "eval_data")
    )
  } else {
    ensemble_metrics <- cbind(as.data.frame(t(fit_ens_eval)), what = "full", origin = "fit_data")
  }

  # Get list of metrics
  metrics_list <- lapply(good_models, function(id) {
    model_fits[[id]]$cv_metrics
  })

  # Convert the list of data frames into a 3-dimensional array
  metrics_array <- array(unlist(metrics_list), dim = c(nrow(metrics_list[[1]]), ncol(metrics_list[[1]]), length(metrics_list)))

  # Calculate the mean across the third dimension (which represents the data frames)
  metrics_mean_df <- apply(metrics_array, c(1, 2), mean, na.rm = T)
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
    metrics_list <- lapply(good_models, function(id) {
      model_fits[[id]]$eval_metrics
    })
    metrics_list <- do.call("rbind", metrics_list)

    # Calculate the mean across
    eval_metrics_mean_df <- apply(metrics_list, 2, mean, na.rm = T)
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

  arrow::write_parquet(
    ensemble_metrics,
    paste0(
      metric_out,
      "taxonid=", species, "_model=", outacro,
      "_method=ensemble_what=cvmetrics.parquet"
    )
  )

  # Ensemble response curves
  to_ensemble <- list.files(metric_out,
                                    pattern = "what=resp", full.names = T)
  to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]
  to_ensemble <- lapply(to_ensemble, arrow::read_parquet)
  
  init_data <- to_ensemble[[1]]
  
  for (k in seq_len(length(to_ensemble))) {
    tdf <- data.frame(to_ensemble[[2]][,c("response")])
    colnames(tdf) <- paste0("response_", k)
    init_data <- cbind(init_data, tdf)
  }
  
  ens_resp <- apply(init_data[,startsWith(colnames(init_data), "resp")], 1, "median")
  ens_resp_sd <- apply(init_data[,startsWith(colnames(init_data), "resp")], 1, "sd")
  
  if ("in_range" %in% colnames(init_data)) {
    ens_resp_curves <- cbind(init_data[,c("variable", "base", "in_range")], response = ens_resp, response_sd = ens_resp_sd)
  } else {
    ens_resp_curves <- cbind(init_data[,c("variable", "base")], response = ens_resp, response_sd = ens_resp_sd)
  }
  
  outfile <- paste0(metric_out,
                    "taxonid=", species, "_model=", outacro, "_method=ensemble",
                    "_what=respcurves.parquet")
  
  arrow::write_parquet(ens_resp_curves, outfile)

  return(list(
    preds = model_predictions,
    logs = model_log
  ))
}



# Save binarization info -----
.cm_save_bin_info <- function(model_predictions, sp_data, algorithms, good_models,
                              metric_out, species, outacro) {
  # Get thresholds
  model_predictions <- lapply(model_predictions, function(x) {x/100})
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
  
  arrow::write_parquet(thresh, paste0(metric_out,
                                "taxonid=", species, "_model=", outacro,
                                "_what=thresholds.parquet"))

  return(thresh_p10_mtp)
}


# Save masks ------
.cm_save_masks <- function(ecoreg_occ, ecoreg_sel, multi_mod, 
                           env, model_predictions, sp_data,
                           pred_out, species, outacro, coord_names,
                           max_depth) {
  # Get base raster
  base_layer <- rast("data/env/terrain/bathymetry_mean.tif")
  bath_layer <- base_layer
  base_layer[!is.na(base_layer)] <- 1

  # Get mask for study area species occurrence
  ecoreg_occ_mask <- terra::mask(base_layer, ecoreg_occ)

  # Get mask for study area ecoregions
  ecoreg_mask <- terra::mask(base_layer, ecoreg_sel)

  # Get area used for fitting
  fit_mask <- terra::extend(env$layers[[1]], model_predictions[[1]])
  fit_mask[!is.na(fit_mask)] <- 1

  # Fitting area + max_depth
  bath_layer <- terra::classify(bath_layer,
                                matrix(c(-999999, max_depth, NA), ncol = 3),
                                right = FALSE)
  max_depth_mask <- terra::mask(fit_mask, bath_layer)

  # Convex hull mask
  conv_hull <- terra::convHull(terra::vect(sp_data$coord_training[sp_data$training$presence == 1, ],
    geom = coord_names, crs = "EPSG:4326"
  ))
  conv_hull_mask <- mask(base_layer, conv_hull)

  minb_circle <- terra::minCircle(terra::vect(sp_data$coord_training[sp_data$training$presence == 1, ],
    geom = coord_names, crs = "EPSG:4326"
  ))
  minb_circle_mask <- mask(base_layer, minb_circle)

  # Buffer mask
  buff_pts <- terra::buffer(
    terra::vect(sp_data$coord_training[sp_data$training$presence == 1, ],
      geom = coord_names, crs = "EPSG:4326"
    ),
    width = 100000
  )
  buff_pts_mask <- terra::mask(base_layer, buff_pts)

  masks <- c(ecoreg_occ_mask, ecoreg_mask, fit_mask, max_depth_mask,
             conv_hull_mask, minb_circle_mask, buff_pts_mask)

  base_layer[!is.na(base_layer)] <- 0

  masks <- terra::mask(base_layer, masks, updatevalue = 1, inverse = T)
  names(masks) <- c(
    "native_ecoregions", "fit_ecoregions", "fit_region", "fit_region_max_depth",
    "convex_hull", "minbounding_circle", "buffer100m"
  )

  masks <- terra::as.int(masks)

  outmask <- paste0(
    pred_out,
    "taxonid=", species, "_model=", outacro, "_mask.tif"
  )
  terra::writeRaster(masks, outmask, overwrite = T, datatype = "INT1U")
  mask_opt <- cogeo_optim(outmask)
  if (file.exists(paste0(outmask, ".aux.json"))) fs::file_delete(paste0(outmask, ".aux.json"))
  return(invisible(NULL))
}


# Post evaluation ----
# Check post-eval number
.cm_check_posteval_n <- function(sp_data, model_log) {
  if (sum(sp_data$training$presence) > 5000) {
    sp_data_post <- sp_data
    which_pres <- which(sp_data$training$presence == 1)
    which_bckg <- which(sp_data$training$presence == 0)
    which_pres <- sample(which_pres, size = 5000)
    sp_data_post$training <- sp_data$training[c(which_pres, which_bckg), ]
    sp_data_post$coord_training <- sp_data$coord_training[c(which_pres, which_bckg), ]
    model_log$other_details <- c(
      model_log$other_details,
      "post evaluation based on 5000 occurrence records"
    )
  } else {
    sp_data_post <- sp_data
  }
  return(list(
    dat = sp_data_post,
    logs = model_log
  ))
}

# Post evaluation temperature
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
  
  binary_maps <- lapply(seq_len(length(model_predictions)), function(x){
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
  
  for (gm in seq_len(length(model_predictions))) {
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
    n_vars <- 2
  }
  
  vars <- vars$variable[seq_len(n_vars)]
  
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
    
    for (pt in seq_len(nrow(var_comb))) {
      
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
      paste0("plot_list[[", seq_len(length(plot_list)), "]]", collapse = "+"), "+ patchwork::plot_layout(guides = 'collect')"
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
    bckg_points <- bckg_points[sample(seq_len(nrow(bckg_points)), 50000),]
  } else if (nrow(bckg_points) < 10000) {
    warning("Low number of background points for niche ecospat. Sample retrieved.")
    env_s <- terra::subset(env_layers, colnames(sdm_data$training)[-1])
    env_s <- as.data.frame(env_s)
    env_s <- na.omit(env_s)
    bckg_points <- env_s[sample(seq_len(nrow(env_s)), 10000),]
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
  
  for (i in seq_len(iterations)) {
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


.cm_posteval_all <- function(post_eval, model_predictions, sp_data_post,
                             thresh_p10_mtp, algorithms, hab_depth,
                             good_models, model_log, metric_out,
                             species, outacro,
                             model_varimport, env, verb_1) {
  
  if ("sst" %in% post_eval) {
    if (verb_1) cli::cli_inform(c(">" = "Performing SST post-evaluation"))
    sst_pe_try <- try(
      .cm_posteval_sst(
        model_predictions, sp_data_post,
        thresh_p10_mtp, algorithms, hab_depth,
        good_models, model_log
      )
    )
    if (!inherits(sst_pe_try, "try-error")) {
      model_log <- sst_pe_try
    }
  }
  if ("niche" %in% post_eval) {
    if (verb_1) cli::cli_inform(c(">" = "Performing niche post-evaluation (ecospat)"))
    if (length(model_predictions) == length(good_models)) {
      new_names <- algorithms[good_models]
    } else {
      new_names <- c(algorithms[good_models], "ensemble")
    }

    niche_eval <- try(lapply(seq_len(length(model_predictions)), function(id) {
      result <- .cm_posteval_nicheequiv(sp_data_post, model_predictions[[id]], env$layers,
        iterations = 5, plot_example = F
      )
      result$model <- new_names[id]
      result
    }), silent = T)
    # We add this as the function is failing some times on the first run
    if (inherits(niche_eval, "try-error")) {
      niche_eval <- try(lapply(seq_len(length(model_predictions)), function(id) {
        result <- .cm_posteval_nicheequiv(sp_data_post, model_predictions[[id]], env$layers,
          iterations = 5, plot_example = F,
          extend.extent = c(-5, 5, -5, 5)
        )
        result$model <- new_names[id]
        result
      }), silent = T)
    }
    if (!inherits(niche_eval, "try-error")) {
      niche_eval <- do.call("rbind", niche_eval)
      arrow::write_parquet(niche_eval, paste0(
        metric_out,
        "taxonid=", species, "_model=", outacro,
        "_what=posteval_niche.parquet"
      ))
      niche_eval <- niche_eval %>%
        group_by(model) %>%
        summarise(across(1:(ncol(.) - 1), function(x) mean(x, na.rm = T)))
      model_log$model_posteval$niche <- as.data.frame(niche_eval[, 1:3])
      # rm(niche_eval)
    } else {
      warning("ecospat niche post-evaluation failed with status\n", niche_eval)
    }
  }
  if ("hyper" %in% post_eval) {
    if (verb_1) cli::cli_inform(c(">" = "Performing niche post-evaluation (hypervolume)"))
    niche_eval <- try(lapply(seq_len(length(model_predictions)), .cm_posteval_hypervolume,
      var_imp = model_varimport,
      model_predictions = model_predictions,
      env_layers = env$layers,
      sdm_data = sp_data_post, return_plot = F
    ), silent = T)
    if (sink.number() > 0) sink()
    if (!inherits(niche_eval, "try-error")) {
      niche_eval_status <- lapply(niche_eval, function(x) {
        r <- x$overlap
        names(r) <- paste0("hyperniche_", names(r))
        data.frame(as.list(r))
      })
      niche_eval_status <- lapply(seq_len(length(niche_eval_status)), function(x) {
        niche_eval_status[[x]]$model <- c(algorithms[good_models], "ensemble")[x]
        niche_eval_status[[x]]
      })
      niche_eval_status <- do.call("rbind", niche_eval_status)

      # Add info to the log
      model_log$model_posteval$hyperniche <- niche_eval_status

      # Save results
      niche_eval_table <- lapply(seq_len(length(model_predictions)), function(x) {
        salg <- c(algorithms[good_models], "ensemble")[x]
        r <- niche_eval[[x]]$occupancy
        r$model <- salg
        r
      })

      niche_eval_table <- dplyr::bind_rows(niche_eval_table)
      niche_eval_table <- dplyr::relocate(niche_eval_table, "model", "occupancy", .before = "type")

      arrow::write_parquet(niche_eval_table, paste0(
        metric_out,
        "taxonid=", species, "_model=", outacro,
        "_what=posteval_hyperniche.parquet"
      ))
      rm(niche_eval, niche_eval_table)
    } else {
      warning("hyperniche post-evaluation failed with status\n", niche_eval)
    }
  }
  return(model_log)
}


# Save metrics -----
.cm_save_metrics <- function(model_fits, metric_out, species,
                             outacro, good_models) {
  lapply(good_models, function(id) {
    model_name <- model_fits[[id]]$name

    outfile <- paste0(
      metric_out,
      "taxonid=", species, "_model=", outacro, "_method=", model_name
    )

    arrow::write_parquet(model_fits[[id]]$cv_metrics, paste0(outfile, "_what=cvmetrics.parquet"))

    other_metrics <- list(
      model_fits[[id]]$full_metrics,
      NULL
    )
    if (!is.null(model_fits[[id]]$eval_metrics)) {
      other_metrics[[2]] <- model_fits[[id]]$eval_metrics
    }
    names(other_metrics) <- c("full", "eval")
    other_metrics <- dplyr::bind_rows(other_metrics, .id = "what")

    arrow::write_parquet(other_metrics, paste0(outfile, "_what=fullmetrics.parquet"))
    return(invisible(NULL))
  })

  return(invisible(NULL))
}


# Save models
.cm_save_models <- function(model_fits, outfolder, species, outacro, good_models) {
  lapply(good_models, function(id) {
    model_name <- model_fits[[id]]$name

    reduced_model <- try(.reduce_model_file(model_fits[[id]]))

    if (!inherits(reduced_model, "try-error")) {
      model_fits[[id]] <- reduced_model
    }

    outpath <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
    if (!dir.exists(outpath)) fs::dir_create(outpath)
    outfile <- paste0(
      outpath,
      "taxonid=", species, "_model=", outacro, "_method=", model_name,
      "_what=model.rds"
    )

    saveRDS(model_fits[[id]], file = outfile)

    return(invisible(NULL))
  })
  return(invisible(NULL))
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
