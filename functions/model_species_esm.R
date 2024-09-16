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
  
  if (algorithms[1] != "esm") {
    algorithms <- "esm"
  }
  
  if (verb_1) cli::cli_alert_info("Starting model for species {species}")
  
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
    model_log$hab_depth <- hab_depth
    
    # Load environmental data
    env <- get_envofgroup(group, depth = hab_depth, load_all = T,
                          conf_file = "sdm_conf.yml", verbose = verb_2)
    
    env <- .cm_check_coastal(species_data, env, coord_names, verb_2)
    
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
    ecoreg_unique <- ecoregions$Realm[is.related(ecoregions,
                                                 vect(bind_rows(fit_pts, eval_pts),
                                                      geom = coord_names),
                                                 "intersects")]

    model_log$model_details$ecoregions <- ecoreg_unique
    ecoreg_occ <- ecoregions[ecoregions$Realm %in% ecoreg_unique,]

    # Apply a buffer to ensure that all areas are covered
    sf::sf_use_s2(FALSE)
    ecoreg_occ_buff <- suppressMessages(
      suppressWarnings(vect(sf::st_buffer(sf::st_as_sf(vect(bind_rows(fit_pts, eval_pts),
                                                            geom = coord_names)), 0.2)))
    )

    adj_ecoreg <- ecoregions$Realm[is.related(ecoregions, ecoreg_occ_buff,
                                              "intersects")]

    # Mask areas
    ecoreg_sel <- ecoregions[ecoregions$Realm %in% unique(
      c(ecoreg_unique, adj_ecoreg)
    ),]
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
    
    bath_pts <- terra::extract(bath, bind_rows(fit_pts, eval_pts))
    
    # Limit by depth if TRUE
    if (limit_by_depth) {
      if (verb_1) cli::cli_alert_info("Limiting by depth")
      
      bath_range <- range(bath_pts[,2])
      bath_range[1] <- bath_range[1] - depth_buffer
      bath_range[2] <- ifelse((bath_range[2] + depth_buffer) > 0, 
                              0, bath_range[2] + depth_buffer)
      
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
    
    
    # Check if species is from shallow/coastal areas and remove distcoast/bathymetry
    if (min(bath_pts[,2], na.rm = T) >= -100) {
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
    
    
    # Assess SAC
    fit_pts_sac <- try(obissdm::outqc_sac_mantel(fit_pts, 
                                                 env_layers = terra::subset(env$layers, env$hypothesis[[1]]),
                                                 plot_result = FALSE,
                                                 verbose = verb_2))
    
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
        model_log$other_details <- c(model_log$other_details,
                                     paste("fit points sampled to", est_tsize))
      } else {
        quad_samp <- est_tsize
      }
    }
    quad_n <- ifelse(env_size_t < quad_samp, round(env_size_t), quad_samp)
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
        
      treg <- obissdm::.get_time(treg, "Sample bias assessment")
    }
    
    # PART 4: MODEL SELECTION ----
    model_log$model_details$hypothesis_tested <- env$hypothesis
    model_log$model_details$best_hypothesis <- names(env$hypothesis)
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
        model_predictions <- lapply(1:length(model_fits), function(id) {
          
          model_name <- paste0("esm_part_", sprintf("%02d", c(1:length(model_fits))[id]))
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
          
          for (k in 1:nrow(scenarios)) {
            
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
            #cogeo_optim(outmess)
          }
          if (do_shape) {
            writeRaster(pred_shape, outshape, overwrite = T, datatype = "INT2U")
            #cogeo_optim(outshape)
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
        
        for (un in 1:length(to_ensemble_g)) {
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
          eval_ens <- extract(ensemble_curr, sp_data$coord_eval, ID = F)
          
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
        
        # Convex hull mask
        conv_hull <- terra::convHull(terra::vect(sp_data$coord_training[sp_data$training$presence == 1,],
                                                 geom = coord_names, crs = "EPSG:4326"))
        conv_hull_mask <- mask(base_layer, conv_hull)
        
        minb_circle <- terra::minCircle(terra::vect(sp_data$coord_training[sp_data$training$presence == 1,],
                                                 geom = coord_names, crs = "EPSG:4326"))
        minb_circle_mask <- mask(base_layer, minb_circle)
        
        # Buffer mask
        buff_pts <-  terra::buffer(terra::vect(sp_data$coord_training[sp_data$training$presence == 1,],
                                                   geom = coord_names, crs = "EPSG:4326"),
                                    width = 100000)
        buff_pts_mask <- mask(base_layer, buff_pts)
        
        masks <- c(ecoreg_occ_mask, ecoreg_mask, fit_mask, conv_hull_mask, minb_circle_mask, buff_pts_mask)
        
        base_layer[!is.na(base_layer)] <- 0
        
        masks <- mask(base_layer, masks, updatevalue = 1, inverse = T)
        names(masks) <- c("native_ecoregions", "fit_ecoregions", "fit_region",
                          "convex_hull", "minbounding_circle", "buffer100m")
        
        masks <- as.int(masks)
        
        outmask <- paste0(pred_out, 
                          "taxonid=", species, "_model=", outacro, "_mask.tif")
        terra::writeRaster(masks, outmask, overwrite = T, datatype = "INT1U")
        mask_opt <- cogeo_optim(outmask)
        if (file.exists(paste0(outmask, ".aux.json"))) fs::file_delete(paste0(outmask, ".aux.json"))
        
        treg <- obissdm::.get_time(treg, "Masks")
        
        
        
        # PART 11: EVALUATE MODELS USING OTHER TECHNIQUES ----
        if (verb_1) cli::cli_alert_info("Performing post-evaluation")
        model_predictions <- list(ensemble_curr)
        good_models <- 1
        
        if ("sst" %in% post_eval) {
          if (verb_1) cli::cli_inform(c(">" = "Performing SST post-evaluation"))
          model_log <- .cm_posteval_sst(model_predictions, sp_data,
                                        thresh_p10_mtp, algorithms, hab_depth,
                                        good_models, model_log)
        }
        if ("niche" %in% post_eval) {
          if (verb_1) cli::cli_inform(c(">" = "Performing niche post-evaluation (ecospat)"))
          if (length(model_predictions) == length(good_models)) {
            new_names <- algorithms[good_models]
          } else {
            new_names <- c(algorithms[good_models], "ensemble")
          }

          niche_eval <- try(lapply(1:length(model_predictions), function(id){
            result <- .cm_posteval_nicheequiv(sp_data, model_predictions[[id]], env$layers,
                                              iterations = 5, plot_example = F)
            result$model <- new_names[id]
            result
          }), silent = T)
          # We add this as the function is failing some times on the first run
          if (inherits(niche_eval, "try-error")) {
            niche_eval <- try(lapply(1:length(model_predictions), function(id){
              result <- .cm_posteval_nicheequiv(sp_data, model_predictions[[id]], env$layers,
                                                iterations = 5, plot_example = F,
                                                extend.extent = c(-5, 5, -5, 5))
              result$model <- new_names[id]
              result
            }), silent = T)
          }
          if (!inherits(niche_eval, "try-error")) {
            niche_eval <- do.call("rbind", niche_eval)
            arrow::write_parquet(niche_eval, paste0(metric_out,
                                                    "taxonid=", species, "_model=", outacro,
                                                    "_what=posteval_niche.parquet"))
            niche_eval <- niche_eval %>% group_by(model) %>% summarise(across(1:(ncol(.)-1), function(x) mean(x,na.rm = T)))
            model_log$model_posteval$niche <- as.data.frame(niche_eval[,1:3])
            #rm(niche_eval)
          } else {
            warning("ecospat niche post-evaluation failed with status\n", niche_eval)
          }
        }
        if ("hyper" %in% post_eval) {
          if (verb_1) cli::cli_inform(c(">" = "Performing niche post-evaluation (hypervolume)"))
          niche_eval <- try(lapply(1:length(model_predictions), .cm_posteval_hypervolume, 
                                   var_imp = list(model_varimport),
                                   model_predictions = model_predictions,
                                   env_layers = env$layers,
                                   sdm_data = sp_data, return_plot = F), silent = T)
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
            
            arrow::write_parquet(niche_eval_table, paste0(metric_out,
                                             "taxonid=", species, "_model=", outacro,
                                             "_what=posteval_hyperniche.parquet"))
            rm(niche_eval, niche_eval_table)
          } else {
            warning("hyperniche post-evaluation failed with status\n", niche_eval)
          }
        }
        
        treg <- obissdm::.get_time(treg, "Post-evaluation")
        
        
        # PART 12: SAVE MODELS OBJECTS ----
        if (verb_1) cli::cli_alert_info("Saving models objects")
        
        # Update log with best parameters
        model_log$model_bestparams$esm <- bind_rows(
          lapply(
            1:length(model_fits), function(x) cbind(model_fits[[x]]$parameters, part = x)
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
