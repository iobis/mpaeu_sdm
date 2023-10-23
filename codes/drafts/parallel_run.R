# Load species list ----
all_sp_groups <- read.csv(list.files("data", pattern = "all_sp", full.names = T)[1])
# Include SDM groups based on the configuration file
all_sp_groups <- get_listbygroup(all_sp)


# Create a storr for control ----
# This gives us more control over what was and what was not done
flux_ctrl_all <- storr_rds("sdmrun")

# Load habitat data ----
eco_info_all <- read.csv("data/species_ecoinfo.csv")

# Start log file ----
if (!file.exists("data/log/sdm_models_log.yml")) {
  log_file <- log_start("data/log/sdm_models_log.yml")
} else {
  log_file <- "data/log/_logtemp"
}




run_model_parallel <- function(
  species,
  eco_info,
  sp_groups,
  log_file,
  flux_ctrl_name = "sdmrun",
  outfold = "results/species",
  outacro = "mv1_pr",
  skip_done = TRUE,
  algos = c("randomforest"),
  run_mode = "skip",
  min_n = 30,
  iflow_tryeval = TRUE,
  ifhigheval = "ifhelps",
  help_perc = 0.25,
  clean_sac = TRUE,
  max_sac = 50,
  std_habitat = "depthsurf",
  try_other_blocks = TRUE,
  min_size_block = 0.2,
  max_size_block = 20,
  test_multi_hypo = c("basevars"),
  fut_scen = c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
  fut_periods = c("dec50", "dec100"),
  clip_area = ext(-41, 47, 20, 89),
  mask_area_names = c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),
  limit_by_depth = TRUE,
  depth_buffer = 300,
  depth_layer = "data/env/terrain/bathymetry_mean.tif",
  scale_vars = FALSE
  ){
  
  library(obissdm)
  library(arrow)
  library(terra)
  library(parallel)
  library(storr)
  source("functions/get_block_n.R")
  source("functions/model_utils.R")
  set.seed(2023)
  
  flux_ctrl <- storr_rds(flux_ctrl_name)
  
  if (!is.null(mask_area_names)) {
    mregions = mregions::mr_shp("MarineRegions:iho")
    mask_area = mregions[mregions$name %in% mask_area_names,]
  }
  
  depth_layer <- rast(depth_layer)
  
  all_sp <- sp_groups
  
  run_models <- function(species) {
    
    if (skip_done) {
      if (flux_ctrl$exists(species)) {
        mod_status <- flux_ctrl$get(sp)$status
        proceed <- !grepl("done", mod_status)
      } else {
        proceed <- TRUE
      }
    } else {
      proceed <- TRUE
    }
    
    if (proceed) {
      timings <- obissdm::.get_time()
      
      # Prepare a function to retrieve the most recent files
      get_recent_file <- function(files) {
        dates <- gsub("\\/.*.", "", gsub("^(.*?)date=", "", files))
        dates <- as.Date(dates, "%Y%m%d")
        files[which.max(dates)]
      }
      
      species_file <- list.files(paste0("data/species/key=", species), 
                                 recursive = T, full.names = T)
      species_file <- species_file[grepl("stdpts", species_file)]
      
      # Put a control here to only model if the species data was standardized
      if (length(species_file) > 0) {
        
        log_file <- log_addsp(log_file, species)
        
        # Open and prepare data
        # Get species group
        sp_group <- all_sp$sdm_group[all_sp$key == species]
        
        # Get eco info for the species
        mode_life <- eco_info$mode_life[eco_info$taxonID == species]
        
        if (mode_life == "pelagic" | mode_life == "pelagic_surface") {
          depth_env <- "depthsurf"
        }
        if (mode_life == "pelagic_mean") {
          depth_env <- "depthmean"
        }
        if (mode_life == "benthic" | mode_life == "demersal" | mode_life == "pelagic_bottom") {
          depth_env <- "depthmax"
        }
        if (mode_life == "NOT_FOUND") {
          depth_env <- std_habitat
        }
        
        # Open species data
        lldim <- c("decimalLongitude", "decimalLatitude")
        
        sp_data <- get_recent_file(species_file)
        sp_data <- read_parquet(sp_data)
        
        # Open environmental data
        hypot <- test_multi_hypo[1]
        
        current <- get_envofgroup(group = sp_group,
                                  depth = depth_env,
                                  verbose = FALSE)
        
        # Make adjustments
        if (!is.null(clip_area)) {
          current <- crop(current, clip_area)
          if (!is.null(mask_area)) {
            current <- mask(current, mask_area, inverse = T)
          }
        }
        if (limit_by_depth) {
          if (!is.null(clip_area)) {
            depth_layer <- crop(depth_layer, clip_area)
            if (!is.null(mask_area)) {
              depth_layer <- mask(depth_layer, mask_area, inverse = T)
            }
          }
          max_depth <- terra::minmax(depth_layer)[1]
          sp_depth <- terra::extract(depth_layer, sp_data[,lldim], ID = F)
          
          sp_max_depth <- ceiling(min(sp_depth[,1], na.rm = T) - depth_buffer)
          
          if (sp_max_depth > max_depth) {
            depth_layer[depth_layer < sp_max_depth] <- NA
            
            current <- mask(current, depth_layer)
          } else {
            sp_max_depth <- max_depth
          }
        } else {
          sp_max_depth <- NA
        }
        if (scale_vars) {
          c_m <- global(current, mean, na.rm = T)[,1]
          c_sd <- global(current, sd, na.rm = T)[,1]
          current <- (current - c_m) / c_sd
        }
        
        # Open future scenarions environmental data
        scenarios_grid <- expand.grid(fut_scen,
                                      fut_periods, stringsAsFactors = F)
        
        for (i in 1:nrow(scenarios_grid)) {
          eval(parse(text = paste0(
            scenarios_grid[i, 1], "_", scenarios_grid[i, 2],
            "<- get_envofgroup(group = sp_group,
                                 depth = depth_env,
                                 scenario = '", scenarios_grid[i, 1],"',
                                 period = '", scenarios_grid[i, 2],"',
                                 verbose = FALSE)"
          )))
          
          if (!is.null(clip_area)) {
            eval(parse(text = paste0(
              scenarios_grid[i, 1], "_", scenarios_grid[i, 2],
              "<- crop(", scenarios_grid[i, 1], "_", scenarios_grid[i, 2],", clip_area)"
            )))
            if (!is.null(mask_area)) {
              eval(parse(text = paste0(
                scenarios_grid[i, 1], "_", scenarios_grid[i, 2],
                "<- mask(", scenarios_grid[i, 1], "_",
                scenarios_grid[i, 2],", mask_area, inverse = T)"
              )))
            }
            if (limit_by_depth & sp_max_depth > max_depth) {
              eval(parse(text = paste0(
                scenarios_grid[i, 1], "_", scenarios_grid[i, 2],
                "<- mask(", scenarios_grid[i, 1], "_",
                scenarios_grid[i, 2],", depth_layer)"
              )))
            }
          }
          if (scale_vars) {
            eval(parse(text = paste0(
              scenarios_grid[i, 1], "_", scenarios_grid[i, 2], "<- (",
              scenarios_grid[i, 1], "_", scenarios_grid[i, 2], " - c_m) / c_sd"
            )))
          }
        }
        timings <- obissdm::.get_time(timings, "load_prepare_envdata")
        
        # Separate evaluation data
        if (any(grepl("eval", sp_data$data_type))) {
          eval_pts <- sp_data[sp_data$data_type == "eval_points", lldim]
        } else {
          eval_pts <- NULL
        }
        
        fit_pts <- sp_data[sp_data$data_type == "fit_points", lldim]
        
        # Check if all points are within the study area
        # Note that the data was standardized considering the full globe - as
        # for some species we expect to, in the future, model the full area of
        # distribution. That's why this second check is necessary.
        fit_pts_valid <- terra::extract(current[[1]], fit_pts, ID = F)
        fit_pts <- fit_pts[!is.na(fit_pts_valid),]
        
        if (!is.null(eval_pts)) {
          eval_pts_valid <- terra::extract(current[[1]], eval_pts, ID = F)
          eval_pts <- eval_pts[!is.na(eval_pts_valid),]
        }
        
        
        # Check if evaluation dataset is bigger than fit points
        if (any(grepl("eval", sp_data$data_type)) & ifhigheval != "ignore") {
          eval_size <- nrow(eval_pts)
          fit_size <- nrow(fit_pts)
          if (eval_size > fit_size) {
            # If so, proceed with one of the strategies:
            total_pts <- rbind(eval_pts, fit_pts)
            # As those are already in the format 1 per cell, we only need to
            # check for duplicates
            total_pts <- total_pts[!duplicated(total_pts),]
            if (ifhigheval == "incorporate") {
              fit_pts <- total_pts
              eval_pts <- NULL
              rm(total_pts)
            }
            if (ifhigheval == "ifhelps") {
              if (nrow(total_pts) > (nrow(fit_pts) + (nrow(fit_pts) * help_perc))) {
                fit_pts <- total_pts
                eval_pts <- NULL
                rm(total_pts)
              } else {
                rm(total_pts)
              }
            }
          }
        }
        
        # If evaluation dataset still exists, then sample 'absence' points to be able
        # to calculate all metrics - but those should be ignored, and only PO metrics
        # used in the end!
        if (!is.null(eval_pts)) {
          valid_c <- as.data.frame(current[[1]], xy = T)
          pa_pts <- valid_c[sample(1:nrow(valid_c), nrow(eval_pts)),1:2]
          colnames(pa_pts) <- lldim
          pa_pts$presence <- 0
          eval_pts$presence <- 1
          eval_pts <- dplyr::bind_rows(eval_pts, pa_pts)
        }
        
        
        # If TRUE, deal with spatial autocorrelation
        if (clean_sac) {
          fit_pts_sac <- try(outqc_sac(fit_pts,
                                       current,
                                       autocor_maxdist = max_sac,
                                       plot_result = F,
                                       verbose = F))
        }
        if (class(fit_pts_sac)[1] != "try-error") {
          fit_pts <- fit_pts_sac
        }
        timings <- obissdm::.get_time(timings, "prepare_spdata")
        
        # If there's enough points, proceed
        if (nrow(fit_pts) > min_n) {
          # See if the number of quadrature points is viable
          total_av_pts <- nrow(as.data.frame(current[[1]]))
          n_quad_number <- ifelse(total_av_pts < 150000, total_av_pts, 150000)
          
          # Prepare data
          sp_data <- mp_prepare_data(fit_pts, species_id = species,
                                     eval_data = eval_pts,
                                     env_layers = current,
                                     quad_number = n_quad_number)
          
          # Prepare cross-validation blocks
          # Retrieve block size
          block_sizes <- get_conf(what = "blocksizes")$blocksizes
          sel_size <- block_sizes[[depth_env]][[sp_group]]
          sel_size <- round(sel_size / 111, 1) # convert to degrees, approximate
          if (sel_size < min_size_block) {
            sel_size <- min_size_block
          } else if (sel_size > max_size_block) {
            sel_size <- max_size_block
          }
          
          # Prepare a grid
          # Ensure that the extent encompasses the resolution
          tune_blocks <- "spatial_grid"
          
          xmin_ext <- round(ext(current)[1]-0.1, 1)
          ymax_ext <- round(ext(current)[4]+0.1, 1)
          
          ymin_t <- round(ext(current)[3]-0.1, 1)
          test_ymin <- seq(ymax_ext, ymin_t, by = -sel_size)
          ymin_ext <- ifelse(min(test_ymin) > ymin_t, round((min(test_ymin) - sel_size), 1), min(test_ymin))
          
          xmax_t <- round(ext(current)[2]+0.1, 1)
          test_xmax <- seq(xmin_ext, xmax_t, by = sel_size)
          xmax_ext <- ifelse(max(test_xmax) < xmax_t, round((max(test_xmax) + sel_size), 1), max(test_xmax))
          
          block_list <- list(spatial_grid = rast(ext(xmin_ext, xmax_ext, ymin_ext, ymax_ext), resolution = sel_size))
          
          sp_data_blocked <- mp_prepare_blocks(sp_data, method = "manual", 
                                               block_types = c("spatial_grid", "random"),
                                               manual_shp = block_list,
                                               n_iterate = 300,
                                               verbose = F)
          
          # Check if any 0
          n_blocks <- get_blocks_n(sp_data_blocked, "spatial_grid")
          n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
          final_block_size <- sel_size
          
          if (any(n_blocks) & try_other_blocks) {
            
            try_sizes <- c((sel_size - sel_size*0.1), (sel_size + sel_size*0.1),
                           (sel_size - sel_size*0.25), (sel_size + sel_size*0.25),
                           (sel_size - sel_size*0.5), (sel_size + sel_size*0.5),
                           (sel_size - sel_size*0.75), (sel_size + sel_size*0.75))
            try_sizes <- try_sizes[try_sizes >= min_size_block & try_sizes <= max_size_block]
            ind <- 1
            
            while (any(n_blocks) & ind <= length(try_sizes)) {
              
              xmin_ext <- round(ext(current)[1]-0.1, 1)
              ymax_ext <- round(ext(current)[4]+0.1, 1)
              
              ymin_t <- round(ext(current)[3]-0.1, 1)
              test_ymin <- seq(ymax_ext, ymin_t, by = -try_sizes[ind])
              ymin_ext <- ifelse(min(test_ymin) > ymin_t, round((min(test_ymin) - try_sizes[ind]), 1), min(test_ymin))
              
              xmax_t <- round(ext(current)[2]+0.1, 1)
              test_xmax <- seq(xmin_ext, xmax_t, by = try_sizes[ind])
              xmax_ext <- ifelse(max(test_xmax) < xmax_t, round((max(test_xmax) + try_sizes[ind]), 1), max(test_xmax))
              
              new_grid <- list(spatial_grid = rast(ext(xmin_ext, xmax_ext, ymin_ext, ymax_ext), resolution = try_sizes[ind]))
              
              sp_data_blocked <- mp_prepare_blocks(sp_data, method = "manual", 
                                                   block_types = c("spatial_grid", "random"),
                                                   manual_shp = new_grid,
                                                   n_iterate = 300,
                                                   verbose = F)
              
              # TODO: Improve to ensure at least 2 points
              n_blocks <- get_blocks_n(sp_data_blocked, "spatial_grid")
              n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
              ind <- ind+1
              final_block_size <- try_sizes[ind]
            }
            
          } else if (any(n_blocks)) {
            sp_data_blocked$blocks$folds <- sp_data_blocked$blocks$folds[
              names(sp_data_blocked$blocks$folds) != "spatial_grid"]
            tune_blocks <- "random"
          }
          
          sp_data <- sp_data_blocked
          rm(sp_data_blocked)
          
          pred_function <- function(index, name_save){
            name_save <- paste0(name_save, "_", scenarios_grid[index, 1], "_",
                                scenarios_grid[index, 2], ".tif")
            eval(parse(text = paste0(
              "pred_", scenarios_grid[index, 1], "_", scenarios_grid[index, 2],
              " <- predict(sp_model, ", scenarios_grid[index, 1], "_", scenarios_grid[index, 2], ")"
            )))
            eval(parse(text = paste0(
              "writeRaster(pred_", scenarios_grid[index, 1],
              "_", scenarios_grid[index, 2],", name_save, overwrite = T)"
            )))
            return(invisible(NULL))
          }
          
          timings <- obissdm::.get_time(timings, "prepare_modelobj")
          
          # Run all models (if more than 1)
          base_outf <- paste0(outfold, "/key=", species, "/", outacro)
          
          # Check algos and what to do
          algos_todo <- algos
          if (run_mode == "skip") {
            if (flux_ctrl$exists(species)) {
              done_algos <- flux_ctrl$get(as.character(species))
              done_algos <- unlist(strsplit(done_algos$algos, "\\|"))
              algos_todo <- algos[!done_algos %in% algos]
            }
          }
          
          # For control
          error_log_list <- list()
          
          if ("lasso_iwlr" %in% algos_todo) {
            # Run model
            sp_model <- try(sdm_module_lasso(sp_data, method = "iwlr",
                                             tune_blocks = tune_blocks))
            # Make predictions
            error_log <- ifok(sp_model, {
              # Prepare outfile
              outf_name <- paste0(
                base_outf,
                "/lasso_iwlr/predictions/", outacro, "_lasso_iwlr_key", species, "_",
                hypot
              )
              fs::dir_create(dirname(outf_name))
              # Predict/save current
              curr_pred <- predict(sp_model, current)
              writeRaster(curr_pred, paste0(outf_name, "_current.tif"), overwrite = T)
              # Predict/save future
              lapply(1:nrow(scenarios_grid), pred_function,
                     name_save = outf_name)
              # Add info to log
              log_addmodel(log_file, "lasso_iwlr", sp_model, sp_data,
                           sp_max_depth, tune_blocks, final_block_size)
              # Save model stuff
              outf_name <- gsub("predictions", "model", outf_name)
              fs::dir_create(dirname(outf_name))
              save_model_stuff(sp_model, current, outf_name)
            })
            error_log_list <- c(error_log_list, lasso_iwlr = error_log)
            timings <- obissdm::.get_time(timings, "model_lasso_iwlr")
          }
          
          if ("lasso_naive" %in% algos_todo) {
            # Run model
            sp_model <- try(sdm_module_lasso(sp_data,
                                             tune_blocks = tune_blocks))
            # Make predictions
            error_log <- ifok(sp_model, {
              # Prepare outfile
              outf_name <- paste0(
                base_outf,
                "/lasso_naive/predictions/", outacro, "_lasso_naive_key", species, "_",
                hypot
              )
              fs::dir_create(dirname(outf_name))
              # Predict/save current
              curr_pred <- predict(sp_model, current)
              writeRaster(curr_pred, paste0(outf_name, "_current.tif"), overwrite = T)
              # Predict/save future
              lapply(1:nrow(scenarios_grid), pred_function,
                     name_save = outf_name)
              # Add info to log
              log_addmodel(log_file, "lasso_naive", sp_model, sp_data,
                           sp_max_depth, tune_blocks, final_block_size)
              # Save model stuff
              outf_name <- gsub("predictions", "model", outf_name)
              fs::dir_create(dirname(outf_name))
              save_model_stuff(sp_model, current, outf_name)
            })
            error_log_list <- c(error_log_list, lasso_naive = error_log)
            timings <- obissdm::.get_time(timings, "model_lasso_naive")
          }
          
          if ("maxent" %in% algos_todo) {
            # Run model
            sp_model <- try(sdm_module_maxent(sp_data,
                                              tune_blocks = tune_blocks))
            # Make predictions
            error_log <- ifok(sp_model, {
              # Prepare outfile
              outf_name <- paste0(
                base_outf,
                "/maxent/predictions/", outacro, "_maxent_key", species, "_",
                hypot
              )
              fs::dir_create(dirname(outf_name))
              # Predict/save current
              curr_pred <- predict(sp_model, current)
              writeRaster(curr_pred, paste0(outf_name, "_current.tif"), overwrite = T)
              # Predict/save future
              lapply(1:nrow(scenarios_grid), pred_function,
                     name_save = outf_name)
              # Add info to log
              log_addmodel(log_file, "maxent", sp_model, sp_data,
                           sp_max_depth, tune_blocks, final_block_size)
              # Save model stuff
              outf_name <- gsub("predictions", "model", outf_name)
              fs::dir_create(dirname(outf_name))
              save_model_stuff(sp_model, current, outf_name)
            })
            error_log_list <- c(error_log_list, maxent = error_log)
            timings <- obissdm::.get_time(timings, "model_maxent")
          }
          
          if ("brt" %in% algos_todo) {
            # Run model
            sp_model <- try(sdm_module_brt(sp_data,
                                           tune_blocks = tune_blocks))
            # Make predictions
            error_log <- ifok(sp_model, {
              # Prepare outfile
              outf_name <- paste0(
                base_outf,
                "/brt/predictions/", outacro, "_brt_key", species, "_",
                hypot
              )
              fs::dir_create(dirname(outf_name))
              # Predict/save current
              curr_pred <- predict(sp_model, current)
              writeRaster(curr_pred, paste0(outf_name, "_current.tif"), overwrite = T)
              # Predict/save future
              lapply(1:nrow(scenarios_grid), pred_function,
                     name_save = outf_name)
              # Add info to log
              log_addmodel(log_file, "brt", sp_model, sp_data,
                           sp_max_depth, tune_blocks, final_block_size)
              # Save model stuff
              outf_name <- gsub("predictions", "model", outf_name)
              fs::dir_create(dirname(outf_name))
              save_model_stuff(sp_model, current, outf_name)
            })
            error_log_list <- c(error_log_list, brt = error_log)
            timings <- obissdm::.get_time(timings, "model_brt")
          }
          
          if ("randomforest" %in% algos_todo) {
            # Run model
            sp_model <- try(sdm_module_rf(sp_data,
                                          tune_blocks = tune_blocks))
            # Make predictions
            error_log <- ifok(sp_model, {
              # Prepare outfile
              outf_name <- paste0(
                base_outf,
                "/rf/predictions/", outacro, "_rf_key", species, "_",
                hypot
              )
              fs::dir_create(dirname(outf_name))
              # Predict/save current
              curr_pred <- predict(sp_model, current)
              writeRaster(curr_pred, paste0(outf_name, "_current.tif"), overwrite = T)
              # Predict/save future
              lapply(1:nrow(scenarios_grid), pred_function,
                     name_save = outf_name)
              # Add info to log
              log_addmodel(log_file, "rf_classification_ds", sp_model, sp_data,
                           sp_max_depth, tune_blocks, final_block_size)
              # Save model stuff
              outf_name <- gsub("predictions", "model", outf_name)
              fs::dir_create(dirname(outf_name))
              save_model_stuff(sp_model, current, outf_name)
            })
            error_log_list <- c(error_log_list, randomforest = error_log)
            timings <- obissdm::.get_time(timings, "model_rf")
          }
          
          flux_ctrl$set(species, data.frame(
            species = species, status = "done", algos = paste(algos, collapse = "|"),
            timings = paste0(paste0(names(timings), "-", round(timings, 2)), collapse = ";"),
            errorlog = paste0(paste0(names(error_log_list), "-", unlist(error_log_list)), collapse = ";")
          ))
          
        } else{
          flux_ctrl$set(species, data.frame(
            species = species, status = "failed-by-lownpoints", algos = paste(algos, collapse = "|"),
            timings = NA,
            errorlog = NA
          ))
        }
      } else {
        flux_ctrl$set(species, data.frame(
          species = species, status = "failed-by-nodata", algos = paste(algos, collapse = "|"),
          timings = NA,
          errorlog = NA
        ))
      }
      
    }
    
    return(invisible(NULL))
  }
  return(invisible(NULL))
}





### Run in parallel
done_list <- storr::storr_rds("datastd")
done_list <- done_list$list()
#result <- mclapply(done_list, run_models, mc.cores = 5)
cl <- makeCluster(2)
result <- parLapply(cl, done_list[1:2], run_model_parallel,
                    eco_info = eco_info_all, sp_groups = all_sp_groups,
                    log_file = log_file)


