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
library(parallel)
source("functions/load_env.R")
source("functions/get_block_n.R")
set.seed(2023)


# Prepare settings ----
outfolder <- "results/vsp_testing/"
fs::dir_create(outfolder)

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
curr <- crop(curr, exp_starea)

# Remove Red sea
mregions <- mregions::mr_shp("MarineRegions:iho")
mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

curr <- mask(curr, mregions, inverse = T)

# Mask to a depth of up to 1500 m
bath_mask <- curr$bathymetry_mean
bath_mask[bath_mask < -1500] <- NA

curr <- mask(curr, bath_mask)

ssp1 <- mask(crop(ssp1, curr[[1]]), curr[[1]])
ssp5 <- mask(crop(ssp5, curr[[1]]), curr[[1]])

# Get the study area total area (in km2) for the DWPR
ar <- terra::cellSize(curr[[1]], unit = "km")
ar <- mask(ar, curr[[1]])
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
                             env_layers = curr,
                             quad_number = 100000
                             #, pred_quad = pred_quad_pts # if using uniform quadrature
                             )
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", manual_shp = block_list,
                               n_iterate = 300, verbose = F)
  
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
      
      write.csv(data.frame(
        metric = names(model$full_metrics),
        value = model$full_metrics
      ), paste0(outfolder, key, "/", name,
                                           "_fullmetrics.csv"), row.names = F)
      write.csv(data.frame(
        metric = names(model$eval_metrics),
        value = model$eval_metrics
      ), paste0(outfolder, key, "/", name,
                                           "_evalmetrics.csv"), row.names = F)
      for (z in names(model$cv_metrics)) {
        write.csv(model$cv_metrics[[z]], paste0(outfolder,  key, "/", name, "_",
                                                z, "metrics.csv"), row.names = F)
      }
      
      timings <- data.frame(what = c("tuning", "cv", "evaluate final", 
                                     "evaluate dataset", "predict",
                                     "response curve"),
                            times = c(
                              model$timings,
                              model$timings[4] + pred_time, 
                              model$timings[4] + pred_time + respc_time
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
  pred_save(m, paste(outacro, "lasso_naive", sep = "_"))
  m <- try(sdm_module_lasso(sp_data, blocks_all = T, method = "iwlr"))
  pred_save(m, paste(outacro, "lasso_iwlr", sep = "_"))
  m <- try(sdm_module_lasso(sp_data, blocks_all = T, method = "dwpr", 
                            total_area = dwpr_val))
  pred_save(m, paste(outacro, "lasso_dpwr", sep = "_"))
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
  if (all(file.exists(paste0(outfolder, key, "/", model_names, "_current.tif")))) {
    
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
    
  }
  
  ### Include test for pseudo-absence case (simple GLM, BRT and Maxent)
  curr_pa <- curr[[1]]
  curr_pa[terra::cellFromXY(curr_pa, as.matrix(train_data[,1:2]))] <- NA
  
  pred_quad_pts <- terra::spatSample(curr, nrow(train_data), na.rm = T, xy = T)
  pred_quad_pts <- cbind(decimalLongitude = pred_quad_pts[,1],
                         decimalLatitude = pred_quad_pts[,2],
                         presence = 0,
                         pred_quad_pts[,3:ncol(pred_quad_pts)])
  
  sp_data <- mp_prepare_data(train_data, eval_data, species_id = paste0("species", key),
                             env_layers = curr,
                             pred_quad = pred_quad_pts # if using uniform quadrature
  )
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", manual_shp = block_list,
                               n_iterate = 300, verbose = F)
  
  # Test if any of the lat blocks have a 0 in test data. If TRUE, remove that block method
  n_blocks <- get_blocks_n(sp_data, "spatial_lat")
  n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
  if (any(n_blocks)) {
    sp_data$blocks$folds <- sp_data$blocks$folds[
      names(sp_data$blocks$folds) != "spatial_lat"]
  }
  
  
  if (verbose) cat("Running maxent (pseudo-absence) \n")
  m <- try(sdm_module_maxent(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "maxnet_pseudo", sep = "_"))
  if (verbose) cat("Running brt_normal (pseudo-absence) \n")
  m <- try(sdm_module_brt(sp_data, weight_resp = F, blocks_all = T))
  pred_save(m, paste(outacro, "brt_normal_pseudo", sep = "_"))
  if (verbose) cat("Running glm_normal (pseudo-absence) \n")
  m <- try(sdm_module_glm(sp_data, method = "normal",
                          weight_resp = F, blocks_all = T))
  pred_save(m, paste(outacro, "glm_normal_pseudo", sep = "_"))
  
  
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
             curr = curr, ssp1 = ssp1, ssp5 = ssp5,
             block_list = block_list,
             bias = params$with_bias[pa],
             nsamp = params$samp_type[pa],
             outacro = save_acro, verbose = T,
             mc.cores = 8)
  }
}
