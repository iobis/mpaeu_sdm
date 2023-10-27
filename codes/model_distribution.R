############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### Model the distribution of species ######################



# Load packages ----
library(obissdm)
library(arrow)
library(terra)
library(parallel)
library(storr)
source("functions/get_block_n.R")
source("functions/model_utils.R")
set.seed(2023)



# Settings ----
## General settings ----
# Output folder
outfold <- "results/species"

# Output acronym
# Output will follow the convention:
# 'outfold/key=species/outacro/algoname' >> subfolders: model, predictions
# Each file will be named as outacro_algo_key{sp_key}_{hypothesis}_{what}.*
outacro <- "mv1_pr"

## Models settings ----
# Should "done", "failed" and "partial done" species be skiped?
# This uses the flux_ctrl "storr" - it will look if the target species
# was already tried and had status equal skip
skip_done <- TRUE
skip_failed <- TRUE
skip_partial <- TRUE

# Which algorithms to run
# Availables are: lasso_naive, lasso_iwlr, maxent, brt and randomforest
algos <- c("lasso_naive", "randomforest")

# Should models (algo) already done be: 
# 1) 'skip', 2) 'redo' (overwrite)
run_mode <- "skip"

# Minimum number of points to be considered for running
min_n <- 30


## Data behavior (species/env) ----
# If number is below 30, should it attempt to use evaluation data (if exists)
iflow_tryeval <- TRUE

# In case the number of eval points is higher than the # of training points:
# 1) 'ignore' (keep running), 2) 'ifhelps' (see if including it increase the
# number of training points by help_perc %), 3) 'incorporate' (ignore eval and
# incorporate as training)
ifhigheval <- "ifhelps"
help_perc <- 0.25

# Remove close points using spatial autocorrelation?
clean_sac <- TRUE

# Maximum distance to be allowed to remove points (in kilometers)
max_sac <- 50

# If no information on the habitat is available, what depth to use?
std_habitat <- "depthsurf"


## Cross-validation settings ----
# If for some reason it's not possible to have at least one occ record in each
# block, try to reduce/increase the size of the blocks.
try_other_blocks <- TRUE
# Minimum and maximum size of blocks (in degrees)
min_size_block <- 0.2
max_size_block <- 20


## Environmental settings ----
# If there are multiple variable hypothesis available, whether to test multiple ones
# If a single value is supplied, then no multi hypothesis testing is done
# Else, the writen hypothesis are tested.
test_multi_hypo <- c("basevars")

# Prediction future scenarios
fut_scen <- c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
fut_periods <- c("dec50", "dec100")

# Clip for study area? Set NULL to skip
clip_area <- c(-41, 47, 20, 89)

# Apply mask? Set mask_area = NULL to skip
# Names are those from mregions
mask_area_names <- c("Red Sea", "Gulf of Aqaba", "Gulf of Suez")

# Limit by species depth? If TRUE, it will mask the layer to include only areas
# from depth 0 up to the maximum depth where the species was registered + depth_buffer
# A depth file should be loaded
limit_by_depth <- TRUE
depth_buffer <- 300
depth_layer_path <- "data/env/terrain/bathymetry_mean.tif" # Path to depth layer

# Scale environmental variables?
scale_vars <- FALSE

## Other settings ----
# Other configuration (variables, groups) is available on the following file
# rstudioapi::documentOpen("sdm_conf.yml")



# Load species list ----
all_sp <- read.csv(list.files("data", pattern = "all_sp", full.names = T)[1])
# Include SDM groups based on the configuration file
all_sp <- get_listbygroup(all_sp)



# Create a storr for control ----
# This gives us more control over what was and what was not done
flux_ctrl <- storr_rds("sdmrun_c")



# Load habitat data ----
eco_info <- read.csv("data/species_ecoinfo.csv")



# Start log file ----
if (!file.exists("data/log/sdm_models_log.yml")) {
  log_file <- log_start("data/log/sdm_models_log.yml")
} else {
  log_file <- "data/log/_logtemp"
}



# Create a function for modeling ----
# Because we will run in parallel, we first create a function that will
# perform all steps
run_models <- function(species) {
  
  species <- as.character(species)
  
  if (skip_done | skip_partial | skip_failed) {
    if (flux_ctrl$exists(species)) {
      mod_status <- flux_ctrl$get(species)$status
      if (skip_done) {
        proceed <- !grepl("done", mod_status)
      }
      if (proceed & skip_partial) {
        proceed <- !grepl("partial", mod_status)
      }
      if (proceed & skip_failed) {
        proceed <- !grepl("failed", mod_status)
      }
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
      # Just check to avoid errors
      if (!species %in% eco_info$taxonID) {
        stop("No mode of life information found for species ", species, "\n")
      }
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
        clip_area <- ext(clip_area)
      }
      if (!is.null(mask_area_names)) {
        mregions <- mregions::mr_shp("MarineRegions:iho")
        mask_area <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]
      }
      
      if (!is.null(clip_area)) {
        current <- crop(current, clip_area)
        if (!is.null(mask_area)) {
          current <- mask(current, mask_area, inverse = T)
        }
      }
      if (limit_by_depth) {
        
        depth_layer <- rast(depth_layer_path)
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
        if (nrow(eval_pts) < 5) { # Check size again, after all
          eval_pts <- NULL
        } else {
          valid_c <- as.data.frame(current, xy = T, na.rm = T)
          pa_pts <- valid_c[sample(1:nrow(valid_c), nrow(eval_pts)),1:2]
          colnames(pa_pts) <- lldim
          pa_pts$presence <- 0
          eval_pts$presence <- 1
          eval_pts <- dplyr::bind_rows(eval_pts, pa_pts)
        }
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
                                   quad_number = n_quad_number,
                                   verbose = FALSE)
        
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
        
        flux_status <- ifelse(any(!is.na(error_log_list)), "partial", "done")
        
        flux_ctrl$set(species, data.frame(
          species = species, status = flux_status, algos = paste(algos, collapse = "|"),
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


# Run models ----
## Temp
done_list <- storr::storr_rds("datastd")
done_list <- done_list$list()
##
# Prepare cluster
cl <- makeCluster(7)

# Export objects
clusterEvalQ(cl, library(obissdm))
clusterEvalQ(cl, library(arrow))
clusterEvalQ(cl, library(terra))
clusterEvalQ(cl, library(storr))
clusterEvalQ(cl, source("functions/get_block_n.R"))
clusterEvalQ(cl, source("functions/model_utils.R"))
to_export <- ls()
to_export <- to_export[!to_export %in% c("cl", "done_list", "to_export")]
clusterExport(cl, varlist = to_export)

# Run models
result <- parLapply(cl, done_list, run_models)
#lapply(c(1017389, 1015919, 27392, 273836), run_models)

# Finalize cluster
stopCluster(cl)



# Finalize model ----
# DANGER ZONE
# UNCOMMENT TO PROCEED
# Finalize log file: just do this once you finalized this run and does not intend
# to add new algorithms or re-run species models
# (adding new species is ok)

# log_close("data/log/sdm_models_log.yml")

# Save flux control log
# and destroy it
# flux_ctrl$destroy()
