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
source("functions/load_env.R")
source("functions/components_model_species.R")
source("functions/auxiliary_modelfit.R")
plan(multisession, workers = 4)

# Prepare settings ----
outfolder <- "results/vsp_testing/"
fs::dir_create(outfolder)


fit_parallel <- function(species, variables, verbose = TRUE) {

  terra::terraOptions(memfrac = (0.6/4))

  files <- list.files(file.path("data/virtual_species", paste0("key=", species)),
                      full.names = T)
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

      status <- try(fit_model_vsp(paste0(species, "_", mo, "_", tf), tf_base, tf_pa, variables, outfolder, verbose = verbose), silent = verbose)
    }
  }
  return(invisible(NULL))
}

fit_model_vsp <- function(species, base_file, pa_file, variables, outfolder, verbose) {

  outacro <- "vsp"

  # Algorithms to be used
  algos <- c("maxent", "rf", "xgboost", "brt")
  # Personalized options
  algo_opts <- obissdm::sdm_options()[algos]
  algo_opts$xgboost$gamma <- c(0, 4)
  algo_opts$xgboost$shrinkage <- c(0.1, 0.3)
  algo_opts$brt$method <- "naive"
  algo_opts$brt$n_trees <- seq(2000, 10000, 2000)
  algo_opts$brt$t_depth <- c(1, 3)
  algo_opts$brt$lr <- c(0.01, 0.001)
  algo_opts$brt$weight_resp <- TRUE
  algo_opts$maxent$features <- c("lq", "h")
  algo_opts$maxent$remult <- seq(0.5, 3, by = 0.5)

  monotonic_vars <- c(
      thetao_min = 1,
      thetao_max = -1,
      sws_min = -1,
      so_min = 1,
      o2_min = 1
  )

  # Should areas be masked by the species depth?
  limit_by_depth <- TRUE
  # A buffer to be applied to the depth limitation
  depth_buffer <- 50
  # Assess spatial bias?
  assess_bias <- TRUE
  # Quadrature size
  quad_samp <- 0.01 # 1% of the total number of points
  # Target metric
  tg_metric <- "cbi"
  # Metric threshold
  tg_threshold <- 0.3

  post_eval <- c("sst", "niche")

  if (verbose) {
    verb_1 <- verb_2 <- TRUE
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

  # Define output paths
  metric_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/metrics/")
  pred_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/")
  mod_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/models/")
  fig_out <- paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/figures/")

  if (verb_1) cli::cli_alert_info("Reading data")
  species_data <- arrow::read_parquet(base_file)
  fit_pts <- species_data[,1:2]

  eval_pts <- arrow::read_parquet(pa_file)[,1:3]

  treg <- obissdm::.get_time(treg, "Species data loading")

  model_log$scientificName <- species

  model_log$hab_depth <- "surface"

  env <- load_env(variables)

  env <- list(hypothesis = list(base = names(env)), layers = env)

  # Load Realms shapefile (we call ecoregions for convenience)
  ecoregions <- vect("data/shapefiles/MarineRealms_BO.shp")

  model_log$n_init_points <- nrow(fit_pts)

  treg <- obissdm::.get_time(treg, "Data loading")

  # PART 2: DATA PREPARING ----
  if (verb_1) cli::cli_alert_info("Preparing data")

  # Check which eco-regions are covered by the points
  ecoreg_unique <- ecoregions$Realm[is.related(
      ecoregions,
      vect(bind_rows(fit_pts, eval_pts),
          geom = coord_names
      ),
      "intersects"
  )]

  model_log$model_details$ecoregions <- ecoreg_unique
  ecoreg_occ <- ecoregions[ecoregions$Realm %in% ecoreg_unique, ]

  # Apply a buffer to ensure that all areas are covered
  sf::sf_use_s2(FALSE)
  ecoreg_occ_buff <- suppressMessages(
      suppressWarnings(vect(sf::st_buffer(sf::st_as_sf(vect(bind_rows(fit_pts, eval_pts),
          geom = coord_names
      )), 0.2)))
  )

  adj_ecoreg <- ecoregions$Realm[is.related(
      ecoregions, ecoreg_occ_buff,
      "intersects"
  )]

  # Mask areas
  ecoreg_sel <- ecoregions[ecoregions$Realm %in% unique(
      c(ecoreg_unique, adj_ecoreg)
  ), ]
  model_log$model_details$ecoregions_included <- unique(ecoreg_sel$Realm)

  # # Apply a buffer to ensure that all areas are covered
  # Not necessary anymore, shapefile was adjusted
  # sf::sf_use_s2(FALSE)
  # ecoreg_sel <- suppressMessages(
  #   suppressWarnings(terra::vect(sf::st_buffer(sf::st_as_sf(terra::aggregate(ecoreg_sel)), 0.02)))
  # )
  ecoreg_sel <- terra::crop(ecoreg_sel, terra::ext(-180, 180, -90, 90))

  # Crop by a limited maximum extension
  # TODO: check if limit by extension.
  # ext_pts <- terra::ext(terra::vect(dplyr::bind_rows(fit_pts, eval_pts[eval_pts$presence == 1,]),
  #                     geom = coord_names))
  #
  # ext_pts <- ext_pts + c(5, 5, 5, 5)
  # ext_pts <- terra::intersect(ext_pts, terra::ext(-180, 180, -90, 90))
  #
  # ecoreg_sel <- terra::crop(ecoreg_sel, ext_pts)


  # Load bathymetry layer
  bath <- terra::rast("data/env/terrain/bathymetry_mean.tif")
  bath <- terra::mask(terra::crop(bath, ecoreg_sel), ecoreg_sel)

  bath_pts <- terra::extract(bath, bind_rows(fit_pts[,1:2], eval_pts[,1:2]))

  # Limit by depth if TRUE
  if (limit_by_depth) {
      if (verb_1) cli::cli_alert_info("Limiting by depth")

      bath_range <- range(bath_pts[, 2])
      bath_range[1] <- bath_range[1] - depth_buffer
      bath_range[2] <- ifelse((bath_range[2] + depth_buffer) > 0,
          0, bath_range[2] + depth_buffer
      )

      bath[bath < bath_range[1] | bath > bath_range[2]] <- NA

      if ("coastal" %in% names(env$hypothesis)) {
          europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
          bath <- terra::crop(bath, europe_starea)
          env$layers <- terra::mask(terra::crop(env$layers, ecoreg_sel), bath)
          env$layers <- terra::mask(env$layers, env$layers$wavefetch)
      } else {
          env$layers <- terra::mask(terra::crop(env$layers, ecoreg_sel), bath)
      }

      model_log$model_details$limited_by_depth <- TRUE
      model_log$model_details$depth_buffer <- depth_buffer
  } else {
      if ("coastal" %in% names(env$hypothesis)) {
          europe_starea <- terra::vect("data/shapefiles/mpa_europe_starea_v2.shp")
          ecoreg_sel <- terra::crop(ecoreg_sel, europe_starea)
          env$layers <- terra::mask(env$layers, ecoreg_sel)
      } else {
          env$layers <- terra::mask(env$layers, ecoreg_sel)
      }
  }

  # Assess SAC
  fit_pts_sac <- try(obissdm::outqc_sac_mantel(fit_pts[,1:2],
      env_layers = terra::subset(env$layers, env$hypothesis[[1]]),
      plot_result = FALSE,
      verbose = verb_2
  ))

  if (inherits(fit_pts_sac, "try-error")) {
      fit_pts_sac <- fit_pts
  }

  # Make data object
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
          model_log$other_details <- c(
              model_log$other_details,
              paste("fit points sampled to", est_tsize)
          )
      } else {
          quad_samp <- est_tsize
      }
  }
  quad_n <- ifelse(env_size_t < quad_samp, round(env_size_t), quad_samp)
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

      require(spatstat)

      # Get spatstat statistics and point process for control
      spat_im <- as.im(as.data.frame(aggregate(env$layers[[1]], 10, na.rm = T), xy = T))
      spat_window <- as.owin(spat_im)

      # Define ppp object based on point locations
      sppcords <- sp_data$coord_training[sp_data$training$presence == 1, ]

      if (nrow(sppcords) > 10000) {
          sppcords <- sppcords %>% dplyr::slice_sample(n = 10000)
      }

      spat_ppp <- ppp(
          x = sppcords$decimalLongitude,
          y = sppcords$decimalLatitude,
          window = spat_window
      )

      # calculate envelope around L-hat estimates.
      spat_env_k <- envelope(spat_ppp, Kest, verbose = F)
      spat_env_l <- envelope(spat_ppp, Lest, verbose = F)

      # save information
      if (!dir.exists(metric_out)) fs::dir_create(metric_out)
      outfile <- paste0(
          metric_out, "taxonid=", species, "_model=", outacro,
          "_what=biasmetrics.rds"
      )
      saveRDS(list(k_stat = spat_env_k, l_stat = spat_env_l),
          file = outfile
      )

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

  algo_opts$xgboost$monotone_constraints <- monotonic_vars[match(multi_mod$best_variables, names(monotonic_vars))]
  algo_opts$brt$var_monotone <- algo_opts$xgboost$monotone_constraints

  # add monotonic 2
  algo_opts$xgboost_mono <- algo_opts$xgboost
  algo_opts$brt_mono <- algo_opts$brt
  algo_opts$xgboost$monotone_constraints <- 1
  algo_opts$brt$var_monotone <- 1

  algos <- c(algos, "xgboost_mono", "brt_mono")
  algorithms <- algos

  treg <- obissdm::.get_time(treg, "Model selection")

  # PART 5: FIT MODELS ----
  if (verb_1) cli::cli_alert_info("Starting model fitting")
  if (!is.null(algo_opts)) {
      algo_opts <- algo_opts[algorithms]
  } else {
      algo_opts <- obissdm::sdm_options()
      algo_opts <- algo_opts[algorithms]
  }
  model_fits <- list()
  for (i in 1:length(algorithms)) {
      if (verb_1) cli::cli_alert_info("Fitting {algorithms[i]}")
      model_fits[[i]] <- try(do.call(
          paste0("sdm_module_", gsub(
            "_mono", "", algorithms[i]
          )),
          list(sp_data,
              verbose = verb_2,
              options = algo_opts[[i]]
          )
      ))
      if (inherits(model_fits[[i]], "try-error")) {
          model_log$model_result[[algorithms[i]]] <- "failed"
          model_fits[i] <- list(NULL)
      } else {
          model_log$model_result[[algorithms[i]]] <- "succeeded"
      }
  }

  treg <- obissdm::.get_time(treg, "Model fit")


  if (any(unlist(model_log$model_result) == "succeeded")) {
      if (verb_1) cli::cli_alert_info("Model fitting concluded. Checking if models are good")

      # PART 6: PREDICTION ----

      # Get what models are above threshold
      good_models <- lapply(model_fits, function(model) {
          if (!is.null(model)) {
              cv_res <- model$cv_metrics
              cv_res <- apply(cv_res, 2, mean, na.rm = T)
              the_metric <- cv_res[[tg_metric]]
              if (!is.na(the_metric)) {
                  if (the_metric >= tg_threshold) {
                      if (sum(is.na(model$cv_metrics[[tg_metric]])) >= 3) {
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
          if (verb_1) cli::cli_alert_info("Good models available, predicting.")
          # Predict models
          model_predictions <- lapply(good_models, function(id) {
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

              # Reduced for testings
              scenarios <- data.frame(scenario = c("current", "ssp126", "ssp585"), year = c(NA,"dec100","dec100"))
              # scenarios <- data.frame(scenario = c("current"), year = c(NA))

              for (k in 1:nrow(scenarios)) {
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
                    env_to_pred <- load_env(variables)
                    env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
                  } else if (scenarios$scenario[k] == "ssp126") {
                    env_to_pred <- load_env(variables, scenario = "ssp126")
                    env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
                  } else {
                    env_to_pred <- load_env(variables, scenario = "ssp585")
                    env_to_pred <- terra::mask(env_to_pred, env$layers[[1]])
                  }

                  pred <- predict(model_fits[[id]], env_to_pred)

                  names(pred) <- paste0(scenarios$scenario[k], ifelse(is.null(period), "", paste0("_", period)))

                  pred <- pred * 100
                  pred <- as.int(pred)
                  writeRaster(pred, outpred, overwrite = T, datatype = "INT1U")
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
                      mess_map_t[cellFromXY(mess_map_t, mess_map[, 1:2])] <- mess_map[, 5]
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

          model_log$model_good <- algorithms[good_models]
          model_log$model_good_metric <- tg_metric
          model_log$model_good_threshold <- tg_threshold



          # PART 7: GET RESPONSE CURVES ----
          if (verb_1) cli::cli_alert_info("Getting response curves")
          # Predict models
          model_respcurves <- lapply(good_models, function(id) {
              model_name <- model_fits[[id]]$name

              if (!dir.exists(metric_out)) fs::dir_create(metric_out)
              outfile <- paste0(
                  metric_out,
                  "taxonid=", species, "_model=", outacro, "_method=", model_name,
                  "_what=respcurves.parquet"
              )

              rcurves <- resp_curves(model_fits[[id]],
                  subset(env$layers, multi_mod$best_variables),
                  sdm_data = sp_data
              )

              arrow::write_parquet(rcurves, outfile)

              if (!dir.exists(fig_out)) fs::dir_create(fig_out)
              outfile <- paste0(
                  fig_out,
                  "taxonid=", species, "_model=", outacro, "_method=", model_name,
                  "_what=responsecurves.png"
              )

              p <- plot(rcurves) + ggtitle(paste("AphiaID:", species, "| Model:", outacro, "| Algo:", model_name))
              ggplot2::ggsave(
                  filename = outfile, plot = p,
                  width = 10, height = 8
              )

              return(invisible(NULL))
          })

          treg <- obissdm::.get_time(treg, "Response curves")


          # Variables importance
          model_varimport <- lapply(good_models, function(id) {
              model_name <- model_fits[[id]]$name
              varimport <- obissdm::variable_importance(model_fits[[id]], sdm_data = sp_data)

              outfile <- paste0(
                  metric_out,
                  "taxonid=", species, "_model=", outacro, "_method=", model_name,
                  "_what=varimportance.parquet"
              )

              arrow::write_parquet(varimport, outfile)
              return(varimport)
          })

          if (length(good_models) > 1) {
              model_varimport_means <- lapply(model_varimport, function(x) {
                  x$mean
              })
              model_varimport_sd <- apply(do.call("cbind", model_varimport_means),
                  1, sd,
                  na.rm = T
              )
              model_varimport_sd <- round(model_varimport_sd, 3)
              model_varimport_means <- apply(do.call("cbind", model_varimport_means),
                  1, mean,
                  na.rm = T
              )
              model_varimport_means <- round(model_varimport_means, 3)
              ens_var_imp <- model_varimport[[1]][, 1:2]
              ens_var_imp$mean <- model_varimport_means
              ens_var_imp$sd <- model_varimport_sd
              outfile <- paste0(
                  metric_out,
                  "taxonid=", species, "_model=", outacro, "_method=ensemble",
                  "_what=varimportance.parquet"
              )

              arrow::write_parquet(ens_var_imp, outfile)
              model_varimport <- c(model_varimport, list(ens_var_imp))
          }

          treg <- obissdm::.get_time(treg, "Variable importance")

          # PART 8: ENSEMBLE ----
          if (length(good_models) > 1) {
              if (verb_1) cli::cli_alert_info("Enough number of models, ensembling")

              to_ensemble <- list.files(paste0(outfolder, "/taxonid=", species, "/model=", outacro, "/predictions/"),
                  pattern = "\\.tif$", full.names = T
              )
              to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]

              to_ensemble_g <- stringr::str_extract(to_ensemble, "scen=*.*_")
              to_ensemble_g <- gsub("scen=", "", to_ensemble_g)
              to_ensemble_g <- unique(to_ensemble_g)

              for (un in 1:length(to_ensemble_g)) {
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

                  ensemble_mean <- as.int(ensemble_mean)
                  ensemble_cv <- as.int(ensemble_cv)
                  writeRaster(c(ensemble_mean, ensemble_cv), outens,
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
                  eval_ens <- extract(ensemble_curr, sp_data$coord_eval, ID = F)

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
          }

          # Ensemble response curves
          if (length(good_models) > 1) {
              to_ensemble <- list.files(metric_out,
                  pattern = "what=resp", full.names = T
              )
              to_ensemble <- to_ensemble[grepl(paste0(substr(algorithms, 1, 3)[good_models], collapse = "|"), to_ensemble)]
              to_ensemble <- lapply(to_ensemble, arrow::read_parquet)

              init_data <- to_ensemble[[1]]

              for (k in 1:length(to_ensemble)) {
                  tdf <- data.frame(to_ensemble[[2]][, c("response")])
                  colnames(tdf) <- paste0("response_", k)
                  init_data <- cbind(init_data, tdf)
              }

              ens_resp <- apply(init_data[, startsWith(colnames(init_data), "resp")], 1, "median")
              ens_resp_sd <- apply(init_data[, startsWith(colnames(init_data), "resp")], 1, "sd")

              if ("in_range" %in% colnames(init_data)) {
                  ens_resp_curves <- cbind(init_data[, c("variable", "base", "in_range")], response = ens_resp, response_sd = ens_resp_sd)
              } else {
                  ens_resp_curves <- cbind(init_data[, c("variable", "base")], response = ens_resp, response_sd = ens_resp_sd)
              }

              outfile <- paste0(
                  metric_out,
                  "taxonid=", species, "_model=", outacro, "_method=ensemble",
                  "_what=respcurves.parquet"
              )

              arrow::write_parquet(ens_resp_curves, outfile)
          }

          treg <- obissdm::.get_time(treg, "Ensemble")

          # PART 9: SAVE BINARIZATION INFO ----
          if (verb_1) cli::cli_alert_info("Saving binarization info")
          # Get thresholds
          model_predictions <- lapply(model_predictions, function(x) {
              x / 100
          })
          # P10 threshold
          thresh_p10_mtp <- lapply(model_predictions, function(pred) {
              predv <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1, ], ID = F)[, 1]
              p10 <- ceiling(length(predv) * 0.9)
              data.frame(p10 = rev(sort(predv))[p10], mtp = min(na.omit(predv)))
          })

          # Max Specificity + Sensitivity
          thresh_maxss <- lapply(model_predictions, function(pred) {
              pred_p <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 1, ], ID = F)[, 1]
              pred_a <- raster::extract(pred, sp_data$coord_training[sp_data$training$presence == 0, ], ID = F)[, 1]

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

          arrow::write_parquet(thresh, paste0(
              metric_out,
              "taxonid=", species, "_model=", outacro,
              "_what=thresholds.parquet"
          ))



          # PART 10: CREATE MASKS AND SAVE ----
          if (verb_1) cli::cli_alert_info("Saving masks")

          # Get base raster
          base_layer <- rast("data/env/terrain/bathymetry_mean.tif")
          base_layer[!is.na(base_layer)] <- 1

          # Get mask for study area species occurrence
          ecoreg_occ_mask <- mask(base_layer, ecoreg_occ)

          # Get mask for study area ecoregions
          ecoreg_mask <- mask(base_layer, ecoreg_sel)

          # Get area used for fitting
          if (multi_mod$best_model == "coastal") {
              fit_mask <- terra::extend(env$layers[[1]], base_layer)
          } else {
              fit_mask <- terra::extend(env$layers[[1]], model_predictions[[1]])
          }
          fit_mask[!is.na(fit_mask)] <- 1

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
          buff_pts_mask <- mask(base_layer, buff_pts)

          masks <- c(ecoreg_occ_mask, ecoreg_mask, fit_mask, conv_hull_mask, minb_circle_mask, buff_pts_mask)

          base_layer[!is.na(base_layer)] <- 0

          masks <- mask(base_layer, masks, updatevalue = 1, inverse = T)
          names(masks) <- c(
              "native_ecoregions", "fit_ecoregions", "fit_region",
              "convex_hull", "minbounding_circle", "buffer100m"
          )

          masks <- as.int(masks)

          outmask <- paste0(
              pred_out,
              "taxonid=", species, "_model=", outacro, "_mask.tif"
          )
          terra::writeRaster(masks, outmask, overwrite = T, datatype = "INT1U")
          mask_opt <- cogeo_optim(outmask)
          if (file.exists(paste0(outmask, ".aux.json"))) fs::file_delete(paste0(outmask, ".aux.json"))

          treg <- obissdm::.get_time(treg, "Masks")



          # PART 11: EVALUATE MODELS USING OTHER TECHNIQUES ----
          if (verb_1) cli::cli_alert_info("Performing post-evaluation")

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

          if ("sst" %in% post_eval) {
              if (verb_1) cli::cli_inform(c(">" = "Performing SST post-evaluation"))
              model_log <- .cm_posteval_sst(
                  model_predictions, sp_data_post,
                  thresh_p10_mtp, algorithms, hab_depth,
                  good_models, model_log
              )
          }
          if ("niche" %in% post_eval) {
              if (verb_1) cli::cli_inform(c(">" = "Performing niche post-evaluation (ecospat)"))
              if (length(model_predictions) == length(good_models)) {
                  new_names <- algorithms[good_models]
              } else {
                  new_names <- c(algorithms[good_models], "ensemble")
              }

              niche_eval <- try(lapply(1:length(model_predictions), function(id) {
                  result <- .cm_posteval_nicheequiv(sp_data_post, model_predictions[[id]], env$layers,
                      iterations = 5, plot_example = F
                  )
                  result$model <- new_names[id]
                  result
              }), silent = T)
              # We add this as the function is failing some times on the first run
              if (inherits(niche_eval, "try-error")) {
                  niche_eval <- try(lapply(1:length(model_predictions), function(id) {
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
              niche_eval <- try(lapply(1:length(model_predictions), .cm_posteval_hypervolume,
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
                  niche_eval_status <- lapply(1:length(niche_eval_status), function(x) {
                      niche_eval_status[[x]]$model <- c(algorithms[good_models], "ensemble")[x]
                      niche_eval_status[[x]]
                  })
                  niche_eval_status <- do.call("rbind", niche_eval_status)

                  # Add info to the log
                  model_log$model_posteval$hyperniche <- niche_eval_status

                  # Save results
                  niche_eval_table <- lapply(1:length(model_predictions), function(x) {
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

          treg <- obissdm::.get_time(treg, "Post-evaluation")


          # PART 12: SAVE MODELS OBJECTS ----
          if (verb_1) cli::cli_alert_info("Saving models objects")

          # Save metrics
          save_metrics <- lapply(good_models, function(id) {
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
              #fs::dir_delete(file.path(outfolder, paste0("taxonid=", species)))
          }
      }
  } else {
      if (verb_1) cli::cli_alert_warning("All models failed")
      st_status <- "all_failed"
      if (dir.exists(file.path(outfolder, paste0("taxonid=", species)))) {
          #fs::dir_delete(file.path(outfolder, paste0("taxonid=", species)))
      }
  }
  return(invisible(NULL))
}

env_vars <- c("thetao-max", "thetao-min", "sws-min", "so-min", "o2-min")

result <- future_map(c(1001, 1002, 1003, 1004),
                     fit_parallel, variables = env_vars, verbose = TRUE,
                     .progress = TRUE, .options = furrr_options(seed = T))

species = 1001

fit_model_vsp(paste0(species, "_", mo, "_", tf), tf_base, tf_pa, variables, outfolder, verbose = T)








# List and load environmental layers ----
env_vars <- c("thetao-mean", "so-mean", "po4-mean", "phyc-mean")

curr <- load_env(env_vars, terrain_vars = "bathymetry_mean")
ssp1 <- load_env(env_vars, scenario = "ssp126", terrain_vars = "bathymetry_mean")
ssp5 <- load_env(env_vars, scenario = "ssp585", terrain_vars = "bathymetry_mean")

# Load study area shapefile
starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
exp_starea <-  ext(-41, 47, 20, 89) # ext(starea) +- 5 degrees

# Crop to the expanded area (only the one that will be used for
# sampling the background)
curr_back <- crop(curr, exp_starea)

# Remove Red sea
mregions <- mregions::mr_shp("MarineRegions:iho")
mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

curr_back <- mask(curr_back, mregions, inverse = T)

# Mask to a depth of up to 1500 m
bath_mask <- curr_back$bathymetry_mean
bath_mask[bath_mask < -1500] <- NA

curr_back <- mask(curr_back, bath_mask)

# Get the study area total area (in km2) for the DWPR
ar <- terra::cellSize(curr_back[[1]], unit = "km")
ar <- mask(ar, curr_back[[1]])
ar <- global(ar, sum, na.rm = T)
dwpr_val <- ar$sum


# Define for which species the model will run ----
keys <- 101:104

# Define a grid for the spatial blocks cross-validation ----
block_list <- list(
  spatial_grid = rast(ext(curr), resolution = 5),
  spatial_lat = rast(ext(curr), ncol = 1, nrow = 30)
)

# Get a uniform sample of quadrature points for all species (optional)
# pred_quad_pts <- terra::as.data.frame(curr[[1]], xy = T, na.rm = T)[,1:2]
# colnames(pred_quad_pts) <- c("decimalLongitude", "decimalLatitude")
# pred_quad_pts <- pred_quad_pts[sample(1:nrow(pred_quad_pts), size = 30000, replace = FALSE),]
# pred_quad_pts <- cbind(pred_quad_pts, presence = 0, terra::extract(curr, pred_quad_pts, ID = F))


# Create a function to run everything in parallel ----
test_methods <- function(repli, key,
                         curr_back,
                         curr, ssp1, ssp5,
                         block_list,
                         bias = F, nsamp = "low",
                         outacro = "sdmtest",
                         verbose = F) {
  
  outacro <- paste0(outacro, "_rep", repli)
  
  if (bias) {
    bias <- "bias_"
  } else {
    bias <- NULL
  }
  
  train_data <- read_parquet(paste0("data/virtual_species/key=", key,
                                    "/occurrences_", nsamp, "_", bias,
                                    "rep", repli, ".parquet"))
  eval_data <- read_parquet(paste0("data/virtual_species/key=", key,
                                   "/occurrences_", nsamp, "_", bias,
                                   "rep", repli, "_pa.parquet"))
  
  sp_data <- mp_prepare_data(train_data, eval_data, species_id = paste0("species", key),
                             env_layers = curr_back,
                             quad_number = 150000
                             #, pred_quad = pred_quad_pts # if using uniform quadrature
                             )
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", manual_shp = block_list,
                               n_iterate = 300)
  
  # Test if any of the lat blocks have a 0 in test data. If TRUE, remove that block method
  n_blocks <- get_blocks_n(sp_data, "spatial_lat")
  n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
  if (any(n_blocks)) {
    sp_data$blocks$folds <- sp_data$blocks$folds[
      names(sp_data$blocks$folds) != "spatial_lat"]
  }
  
  
  pred_save <- function(model, name) {
    if (class(model)[1] != "try-error") {
      fs::dir_create(paste0(outfolder, key))
      pred_time <- obissdm:::.get_time()
      pred_curr <- predict(model, curr)
      pred_time <- obissdm:::.get_time(pred_time)
      pred_ssp1 <- predict(model, ssp1)
      pred_ssp5 <- predict(model, ssp5)
      
      respc_time <- obissdm:::.get_time()
      resp_cur <- resp_curves(m, curr)
      write.csv(resp_cur, paste0(outfolder, key, "/", name, "_respcurves.csv"),
                row.names = F)
      respc_time <- obissdm:::.get_time(respc_time)
      
      write.csv(model$full_metrics, paste0(outfolder, key, "/", name,
                                           "_fullmetrics.csv"), row.names = F)
      write.csv(model$eval_metrics, paste0(outfolder, key, "/", name,
                                           "_evalmetrics.csv"), row.names = F)
      for (z in names(model$cv_metrics)) {
        write.csv(model$cv_metrics[[z]], paste0(outfolder,  key, "/", name, "_",
                                                z, "metrics.csv"), row.names = F)
      }
      
      timings <- data.frame(what = c("tuning", "cv", "evaluate final", 
                                     "evaluate dataset", "predict",
                                     "response curve"),
                            times = c(
                              model$timings, pred_time, respc_time
                            ))
      
      write.csv(timings, paste0(outfolder, key, "/", name, "_timings.csv"), 
                row.names = F)
      
      writeRaster(pred_curr, paste0(outfolder, key, "/", name, "_current.tif"),
                  overwrite = T)
      writeRaster(pred_ssp1, paste0(outfolder, key, "/", name, "_ssp1.tif"),
                  overwrite = T)
      writeRaster(pred_ssp5, paste0(outfolder, key, "/", name, "_ssp5.tif"),
                  overwrite = T)
      
    } 
    return(invisible(NULL))
  }
  
  # Maxent
  if (verbose) cat("Running maxent \n")
  m <- try(sdm_module_maxent(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "maxnet", sep = "_"))
  # Lasso
  if (verbose) cat("Running LASSO \n")
  m <- try(sdm_module_lasso(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "lasso", sep = "_"))
  # GLM
  if (verbose) cat("Running glm_normal \n")
  m <- try(sdm_module_glm(sp_data, method = "normal",
                          weight_resp = F, blocks_all = T))
  pred_save(m, paste(outacro, "glm_normal", sep = "_"))
  if (verbose) cat("Running glm_naive \n")
  m <- try(sdm_module_glm(sp_data, method = "normal",
                          weight_resp = T, blocks_all = T))
  pred_save(m, paste(outacro, "glm_naive", sep = "_"))
  if (verbose) cat("Running glm_iwlr \n")
  m <- try(sdm_module_glm(sp_data, method = "iwlr", blocks_all = T))
  pred_save(m, paste(outacro, "glm_iwlr", sep = "_"))
  if (verbose) cat("Running glm_dwpr \n")
  m <- try(sdm_module_glm(sp_data, method = "dwpr",
                          total_area = dwpr_val, blocks_all = T))
  pred_save(m, paste(outacro, "glm_dwpr", sep = "_"))
  # GAM
  # m <- sdm_module_gam(sp_data, method = "normal", weight_resp = F, blocks_all = T)
  # pred_save(m, paste(outacro, "gam_normal", nsamp, bias, sep = "_"))
  # m <- sdm_module_gam(sp_data, method = "normal", weight_resp = T, blocks_all = T)
  # pred_save(m, paste(outacro, "gam_naive", nsamp, bias, sep = "_"))
  # BRT
  if (verbose) cat("Running brt_normal \n")
  m <- try(sdm_module_brt(sp_data, weight_resp = F, blocks_all = T))
  pred_save(m, paste(outacro, "brt_normal", sep = "_"))
  if (verbose) cat("Running brt_naive \n")
  m <- try(sdm_module_brt(sp_data, weight_resp = T, blocks_all = T))
  pred_save(m, paste(outacro, "brt_naive", sep = "_"))
  # RandomForest
  if (verbose) cat("Running rf_downsampled \n")
  m <- try(sdm_module_rf(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "rf_downsampled", sep = "_"))
  if (verbose) cat("Running rf_classification \n")
  m <- try(sdm_module_rf(sp_data, method = "classification",
                         type = "normal", blocks_all = T))
  pred_save(m, paste(outacro, "rf_classification", sep = "_"))
  # LightGBM
  if (verbose) cat("Running lgbm \n")
  m <- try(sdm_module_lgbm(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "lgbm", sep = "_"))
  
  # Make ensemble
  if (verbose) cat("Running ensemble \n")
  model_names <- paste(outacro, c(
    "maxnet", "lasso", "brt_naive", "rf_downsampled"
  ), sep = "_")
  
  # Load predictions for each model
  pred_curr <- rast(paste0(outfolder, key, "/", model_names, "_current.tif"))
  pred_ssp1 <- rast(paste0(outfolder, key, "/", model_names, "_ssp1.tif"))
  pred_ssp5 <- rast(paste0(outfolder, key, "/", model_names, "_ssp5.tif"))
  
  # Scale all to be between 0-1 (remember, the important is the relative value in our models)
  pred_curr <- (pred_curr - global(pred_curr, min, na.rm = T)[,1]) / 
    (global(pred_curr, max, na.rm = T)[,1] - global(pred_curr, min, na.rm = T)[,1])
  pred_ssp1 <- (pred_ssp1 - global(pred_ssp1, min, na.rm = T)[,1]) / 
    (global(pred_ssp1, max, na.rm = T)[,1] - global(pred_ssp1, min, na.rm = T)[,1])
  pred_ssp5 <- (pred_ssp5 - global(pred_ssp5, min, na.rm = T)[,1]) / 
    (global(pred_ssp5, max, na.rm = T)[,1] - global(pred_ssp5, min, na.rm = T)[,1])
  
  pred_curr <- mean(pred_curr)
  pred_ssp1 <- mean(pred_ssp1)
  pred_ssp5 <- mean(pred_ssp5)
  
  writeRaster(pred_curr, paste0(outfolder, key, "/",
                                paste(outacro, "ensemble", sep = "_"),
                                "_current.tif"), overwrite = T)
  writeRaster(pred_ssp1, paste0(outfolder, key, "/",
                                paste(outacro, "ensemble", sep = "_"),
                                "_ssp1.tif"), overwrite = T)
  writeRaster(pred_ssp5, paste0(outfolder, key, "/",
                                paste(outacro, "ensemble", sep = "_"),
                                "_ssp5.tif"), overwrite = T)
  
  # Do the same for the response curve
  resp_cur_files <- paste0(outfolder, key, "/", model_names, "_respcurves.csv")
  resp_curves <- lapply(resp_cur_files, read.csv)
  vars <- unique(resp_curves[[1]][,1])
  for (z in 1:length(vars)) {
    for (k in 1:length(resp_curves)) {
      if (k == 1) {
        to_merge <- resp_curves[[k]][resp_curves[[k]][,1] == vars[z], "response"]
      } else {
        to_merge <- cbind(to_merge,
                          resp_curves[[k]][resp_curves[[k]][,1] == vars[z], "response"])
      }
    }
    merged <- apply(to_merge, 1, mean, na.rm = T)
    new_resp <- data.frame(variable = vars[z],
                      response = merged,
                      base = resp_curves[[1]][resp_curves[[1]][,1] == vars[z], "base"])
    if (z == 1) {
      ens_resp_curves <- new_resp
    } else {
      ens_resp_curves <- rbind(ens_resp_curves, new_resp)
    }
  }
  
  write.csv(ens_resp_curves, paste0(outfolder, key, "/", 
                                    paste(outacro, "ensemble", sep = "_"),
                                    "_respcurves.csv"), row.names = F)
  
  return(invisible(NULL))
}


# Run testing ----
params <- expand.grid(with_bias = c(TRUE, FALSE),
                      samp_type = c("low", "high"))


for (sp in keys) {
  for (pa in 1:nrow(params)) {
    
    save_acro <- paste0(
      "sp", sp, "_",
      ifelse(params$with_bias[pa], "bias", "nobias"),
      "_", params$samp_type[pa]
    )
    
    mclapply(1:10, test_methods, key = sp,
             curr_back = curr_back,
             curr = curr, ssp1 = ssp1, ssp5 = ssp5,
             block_list = block_list,
             bias = params$with_bias[pa],
             nsamp = params$samp_type[pa],
             outacro = save_acro,
             mc.cores = 8)
  }
}
