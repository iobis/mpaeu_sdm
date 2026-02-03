############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### Model the distribution of species ######################

# NOTE: those tests were done with a previous version of the framework, and
# much was updated since then. Those species were modelled with the new method
# and can be found on the S3 bucket.


# Load packages / settings ----
library(obissdm)
library(arrow)
library(terra)
library(parallel)
source("functions/load_env.R")
source("functions/get_block_n.R")
source("functions/get_hab_info.R")
set.seed(2023)

# Settings ----
outfolder <- "results/species_test/"
fs::dir_create(outfolder)


# Load species list ----
species_list <- c(127138, # Hippoglossus hippoglossus
                  127061, # Pagrus auriga
                  127063, # Pagrus pagrus
                  127160, # Solea solea
                  126440) # Pollachius pollachius



# Create the function with all the steps ----
model_sp <- function(species) {
  
  # Get mode of life
  mode_life <- get_hab_info(species)
  
  if (mode_life == "NOT_FOUND" | mode_life == "pelagic" | mode_life == "pelagic_surface") {
    depth_env <- "surf"
  }
  if (mode_life == "pelagic_mean") {
    depth_env <- "mean"
  }
  if (mode_life == "benthic" | mode_life == "demersal" | mode_life == "pelagic_bottom") {
    depth_env <- "max"
  }
  
  # Load list of environmental variables ----
  env_vars <- c("thetao-mean", "so-mean", "po4-mean", "phyc-mean", "sws-max")
  
  curr <- load_env(env_vars, depth = depth_env, terrain_vars = "bathymetry_mean")
  ssp1 <- load_env(env_vars, depth = depth_env, scenario = "ssp126", terrain_vars = "bathymetry_mean")
  ssp5 <- load_env(env_vars, depth = depth_env, scenario = "ssp585", terrain_vars = "bathymetry_mean")
  
  # Load study area shapefile
  # starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
  exp_starea <-  ext(-41, 47, 20, 89) # ext(starea) +- 5 degrees
  
  # Crop to the expanded area (only the one that will be used for
  # sampling the background)
  
  # Remove Red sea
  mregions <- mregions::mr_shp("MarineRegions:iho")
  mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]
  
  curr <- mask(crop(curr, exp_starea), mregions, inverse = T)
  
  # See all available files
  sp_files <- list.files(paste0("data/species/key=", species),
                         recursive = T, full.names = T)
  
  # Get the most recent one
  dates <- stringr::str_extract(sp_files, "date=[:digit:]*[:digit:]")
  dates <- lubridate::as_date(gsub("date=", "", dates), format = "%Y%m%d")
  sel_date <- min(dates)
  
  sp_files <- sp_files[grepl(format(sel_date, "%Y%m%d"), sp_files)]
  
  # Open file
  #pts <- read_parquet()
  
  #### TEMP
  pts <- lapply(sp_files, function(x){
    file <- read_parquet(x)
    file[,c("decimalLongitude", "decimalLatitude")]
  })
  pts <- do.call("rbind", pts)
  
  # get 1 point per cell
  base <- curr[[1]]
  base[] <- NA
  base[cellFromXY(base, as.data.frame(pts))] <- 1
  base <- mask(base, curr[[1]])
  
  pts <- as.data.frame(base, xy = T)[,1:2]
  colnames(pts) <- c("decimalLongitude", "decimalLatitude")
  ####
  
  depth_sp <- terra::extract(curr$bathymetry_mean, pts, ID = F)[,1]
  depth_lim <- ceiling(min(depth_sp) - 500)
  
  # Mask to depth
  bath_mask <- curr$bathymetry_mean
  bath_mask[bath_mask < depth_lim] <- NA
  
  curr <- mask(curr, bath_mask)
  
  ssp1 <- mask(crop(ssp1, exp_starea), curr)
  ssp5 <- mask(crop(ssp5, exp_starea), curr)
  
  # See if there is a testing dataset designates
  # if ("testing" %in% pts$data_type) {
  #   testing <- pts[pts$data_type == "testing", ]
  # } else {
  #   testing <- NULL
  # }
  # 
  # training <- pts[pts$data_type == "training", ]
  
  # Blocks
  block_list <- list(
    spatial_grid = rast(ext(curr), resolution = 5),
    spatial_lat = rast(ext(curr), ncol = 1, nrow = 30)
  )
  
  # Prepare data for models
  sp_data <- mp_prepare_data(pts, species_id = species,
                             env_layers = curr,
                             quad_number = 150000)
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", 
                               manual_shp = block_list,
                               n_iterate = 300)
  
  # Test if any of the lat blocks have a 0 in test data. If TRUE, remove that block method
  n_blocks <- get_blocks_n(sp_data, "spatial_grid")
  n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
  if (any(n_blocks)) {
    sp_data$blocks$folds <- sp_data$blocks$folds[
      names(sp_data$blocks$folds) != "spatial_lat"]
  }
  
  pred_save <- function(model, name) {
    if (class(model)[1] != "try-error") {
      key <- species
      
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
                                     "predict",
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
  
  # Model
  outacro <- paste0("sp", species, "_sreso")
  
  #cat("Running lasso")
  m <- try(sdm_module_lasso(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "lasso_naive", sep = "_"))
  # m <- try(sdm_module_lasso(sp_data, blocks_all = T, method = "iwlr"))
  # pred_save(m, paste(outacro, "lasso_iwlr", sep = "_"))
  m <- try(sdm_module_maxent(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "maxnet", sep = "_"))
  m <- try(sdm_module_brt(sp_data, weight_resp = T, blocks_all = T))
  pred_save(m, paste(outacro, "brt_naive", sep = "_"))
  m <- try(sdm_module_rf(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "rf_downsampled", sep = "_"))
  
  
  
  #### Aggregated test
  
  #### TEMP
  pts <- lapply(sp_files, function(x){
    file <- read_parquet(x)
    file[,c("decimalLongitude", "decimalLatitude")]
  })
  pts <- do.call("rbind", pts)
  
  # get 1 point per cell
  base <- curr[[1]]
  base <- terra::aggregate(base, fact = 2)
  base[] <- NA
  base[cellFromXY(base, as.data.frame(pts))] <- 1
  base <- mask(base, terra::aggregate(curr[[1]], fact = 2))
  
  pts <- as.data.frame(base, xy = T)[,1:2]
  colnames(pts) <- c("decimalLongitude", "decimalLatitude")
  ####
  
  # See if there is a testing dataset designates
  # if ("testing" %in% pts$data_type) {
  #   testing <- pts[pts$data_type == "testing", ]
  # } else {
  #   testing <- NULL
  # }
  # 
  # training <- pts[pts$data_type == "training", ]
  
  # Prepare data for models
  sp_data <- mp_prepare_data(pts, species_id = species,
                             env_layers = curr,
                             quad_number = 150000)
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", 
                               manual_shp = block_list,
                               n_iterate = 300)
  
  # Test if any of the lat blocks have a 0 in test data. If TRUE, remove that block method
  n_blocks <- get_blocks_n(sp_data, "spatial_lat")
  n_blocks <- apply(n_blocks, 2, function(x) ifelse(x == 0, TRUE, FALSE))
  if (any(n_blocks)) {
    sp_data$blocks$folds <- sp_data$blocks$folds[
      names(sp_data$blocks$folds) != "spatial_lat"]
  }
  
  
  # Model
  outacro <- paste0("sp", species, "_areso")
  
  m <- try(sdm_module_lasso(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "lasso_naive", sep = "_"))
  # m <- try(sdm_module_lasso(sp_data, blocks_all = T, method = "iwlr"))
  # pred_save(m, paste(outacro, "lasso_iwlr", sep = "_"))
  m <- try(sdm_module_maxent(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "maxnet", sep = "_"))
  m <- try(sdm_module_brt(sp_data, weight_resp = T, blocks_all = T))
  pred_save(m, paste(outacro, "brt_naive", sep = "_"))
  m <- try(sdm_module_rf(sp_data, blocks_all = T))
  pred_save(m, paste(outacro, "rf_downsampled", sep = "_"))
  
  
  return(invisible(NULL))
  
}


# Run in parallel ----
# mclapply(species_list, model_sp, mc.cores = 5)

model_sp(species_list[1])


# Compare results with IUCN range
iucn_compare <- data.frame(
  species = NA,
  model = paste(c("sreso", "areso"),
                c("lasso_naive", "lasso_iwlr", "maxnet", "brt_naive", "rf_downsampled"),
                sep = "_"),
  jacc = NA
)

for (sp in species_list) {

  iucdir <- tempdir()

  unzip(paste0("data/iucn/iucn_", sp, ".zip"), exdir = iucdir)

  iucn_range <- vect(paste0(iucdir, "/data_0.shp"))

  iucn_compare$species <- sp




  # Get mode of life
  mode_life <- get_hab_info(sp)

  if (mode_life == "NOT_FOUND" | mode_life == "pelagic" | mode_life == "pelagic_surface") {
    depth_env <- "surf"
  }
  if (mode_life == "pelagic_mean") {
    depth_env <- "mean"
  }
  if (mode_life == "benthic" | mode_life == "demersal" | mode_life == "pelagic_bottom") {
    depth_env <- "max"
  }

  # Load list of environmental variables ----
  env_vars <- c("thetao-mean", "so-mean", "po4-mean", "phyc-mean", "sws-max")

  curr <- load_env(env_vars, depth = depth_env, terrain_vars = "bathymetry_mean")

  # Load study area shapefile
  # starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
  exp_starea <-  ext(-41, 47, 20, 89) # ext(starea) +- 5 degrees

  # Crop to the expanded area (only the one that will be used for
  # sampling the background)

  # Remove Red sea
  mregions <- mregions::mr_shp("MarineRegions:iho")
  mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

  curr <- mask(crop(curr, exp_starea), mregions, inverse = T)

  # See all available files
  sp_files <- list.files(paste0("data/species/key=", sp),
                         recursive = T, full.names = T)

  # Get the most recent one
  dates <- stringr::str_extract(sp_files, "date=[:digit:]*[:digit:]")
  dates <- lubridate::as_date(gsub("date=", "", dates), format = "%Y%m%d")
  sel_date <- min(dates)

  sp_files <- sp_files[grepl(format(sel_date, "%Y%m%d"), sp_files)]

  # Open file
  #pts <- read_parquet()

  #### TEMP
  pts <- lapply(sp_files, function(x){
    file <- read_parquet(x)
    file[,c("decimalLongitude", "decimalLatitude")]
  })
  pts <- do.call("rbind", pts)

  # get 1 point per cell
  base <- curr[[1]]
  base[] <- NA
  base[cellFromXY(base, as.data.frame(pts))] <- 1
  base <- mask(base, curr[[1]])

  pts_full <- as.data.frame(base, xy = T)[,1:2]
  colnames(pts_full) <- c("decimalLongitude", "decimalLatitude")
  ####

  depth_sp <- terra::extract(curr$bathymetry_mean, pts_full[,1:2], ID = F)[,1]
  depth_lim <- ceiling(min(depth_sp) - 500)

  # Mask to depth
  bath_mask <- curr$bathymetry_mean
  bath_mask[bath_mask < depth_lim] <- NA

  curr <- mask(curr, bath_mask)

  #### TEMP
  pts <- lapply(sp_files, function(x){
    file <- read_parquet(x)
    file[,c("decimalLongitude", "decimalLatitude")]
  })
  pts <- do.call("rbind", pts)

  # get 1 point per cell
  base <- curr[[1]]
  base <- terra::aggregate(base, fact = 2)
  base[] <- NA
  base[cellFromXY(base, as.data.frame(pts))] <- 1
  base <- mask(base, terra::aggregate(curr[[1]], fact = 2))

  pts_agg <- as.data.frame(base, xy = T)[,1:2]
  colnames(pts_agg) <- c("decimalLongitude", "decimalLatitude")
  ####

  back <- spatSample(curr[[1]], 1000, xy = T, na.rm = T)[,1:2]
  colnames(back) <- c("decimalLongitude", "decimalLatitude")

  pts_agg <- rbind(cbind(pts_agg, presence = 1), cbind(back, presence = 0))
  pts_full <- rbind(cbind(pts_full, presence = 1), cbind(back, presence = 0))


  iucn_rast <- curr[[1]]

  iucn_rast[] <- 1

  iucn_rast <- mask(iucn_rast, crop(iucn_range, iucn_rast), updatevalue = 0)
  iucn_rast <- mask(iucn_rast, curr[[1]])

  for (z in 1:nrow(iucn_compare)) {
    if (file.exists(paste0("results/species_test/", sp, "/sp", sp, "_", iucn_compare$model[z], "_current.tif"))) {
      pred <- rast(paste0("results/species_test/", sp, "/sp", sp, "_", iucn_compare$model[z], "_current.tif"))
      pred <- (pred - global(pred, min, na.rm = T)[,1]) /
        (global(pred, max, na.rm = T)[,1] - global(pred, min, na.rm = T)[,1])
      
      if (grepl("areso", iucn_compare$model[z])) {
        pred_pts <- terra::extract(pred, pts_agg[,1:2], ID = F)[,1]
        obs <- pts_agg$presence
      } else {
        pred_pts <- terra::extract(pred, pts_full[,1:2], ID = F)[,1]
        obs <- pts_full$presence
      }
      
      th <- modEvA::getThreshold(obs = obs, pred = pred_pts,
                                 threshMethod = "MTP", quant = 0.1)
      
      pred_bin <- pred
      pred_bin[pred_bin < th] <- 0
      pred_bin[pred_bin >= th] <- 1
      
      jacc <- function(rast1, rast2) {
        comb <- rast1 + rast2
        inter <- comb
        union <- comb
        inter[] <- ifelse(union[] >= 2, 1, 0)
        union[] <- ifelse(comb[] >= 1, 1, 0)
        
        cinter <- freq(inter)
        cunion <- freq(union)
        
        return((cinter$count[cinter$value == 1] / cunion$count[cunion$value == 1]))
      }
      
      iucn_compare$jacc[z] <- jacc(pred_bin, iucn_rast)
      
      writeRaster(pred_bin,
                  paste0("results/species_test/", sp, "/sp", sp, "_", iucn_compare$model[z], "_binp10_current.tif"),
                  overwrite = T)
    }
  }

  if (sp == species_list[1]) {
    iucn_allsp <- iucn_compare
  } else {
    iucn_allsp <- rbind(iucn_allsp, iucn_compare)
  }

  iucn_compare$species <- iucn_compare$jacc <- NA

}


iucn_res <- iucn_allsp %>%
  filter(!grepl("iwlr", model)) %>%
  mutate(model = gsub("areso", "AGG", model)) %>%
  mutate(model = gsub("sreso", "STD", model)) %>%
  mutate(model = toupper(gsub("_", " ", model))) %>%
  mutate(species = ifelse(species == 127138, "Hippoglossus hippoglossus",
                          ifelse(species == 127061, "Pagrus auriga",
                                 ifelse(species == 127063, "Pagrus pagrus",
                                        "Solea solea"))))

write.csv(iucn_res, "iucn_jacc_res.csv", row.names = F)


# 
# par(mfrow = c(1,3)) 
# iuc <- tempdir()
# unzip("data/iucn/iucn_127063.zip", exdir = iuc)
# iucn <- sf::st_read(paste0(iuc, "/data_0.shp"))
# pred_lasso <- rast("results/species_test/127063/sp127063_sreso_lasso_naive_current.tif")
# pred_rf <- rast("results/species_test/127063/sp127063_sreso_rf_downsampled_current.tif")
# 
# pred_lasso <- as.data.frame(pred_lasso, xy = T)
# colnames(pred_lasso)[3] <- "values"
# 
# pred_rf <- as.data.frame(pred_rf, xy = T)
# colnames(pred_rf)[3] <- "values"
# 
# 
# base <- rnaturalearth::ne_countries(returnclass = "sf")
# 
# 
# (p1 <- ggplot()+
#     geom_sf(data = base, fill = "grey80", color = "grey80") +
#     geom_sf(data = iucn, fill = "#17AB92") +
#     coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
#     geom_point(data = as.data.frame(pts_full), aes(x = decimalLongitude, y = decimalLatitude),
#                size = .5) +
#     #scale_fill_distiller("", direction = 1) +
#     xlab(NULL) + ylab(NULL) +
#     theme_light() +
#     ggtitle(expression(paste(italic("Pagrus pagrus"), " IUCN range map")),
#             "Points are the occurrence records used for modeling.")+
#     theme(panel.border = element_blank(),
#           legend.position = "none",
#           strip.background = element_blank(),
#           strip.text = element_text(color = "grey10"),
#           text = element_text(size = 6)))
# 
# (p2 <- ggplot()+
#   geom_sf(data = base, fill = "grey80", color = "grey80") +
#   geom_raster(data = pred_lasso, aes(x = x, y = y, fill = values)) +
#   coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
#   scale_fill_distiller("", direction = 1) +
#   xlab(NULL) + ylab(NULL) +
#   theme_light() +
#   ggtitle("LASSO") +
#   theme(panel.border = element_blank(),
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(color = "grey10"),
#         text = element_text(size = 6)))
# 
# (p3 <- ggplot()+
#     geom_sf(data = base, fill = "grey80", color = "grey80") +
#     geom_raster(data = pred_rf, aes(x = x, y = y, fill = values)) +
#     coord_sf(xlim = c(-34, 41), ylim = c(24.5, 84.5)) +
#     scale_fill_distiller("", direction = 1) +
#     xlab(NULL) + ylab(NULL) +
#     theme_light() +
#     ggtitle("RF down-sampled") +
#     theme(panel.border = element_blank(),
#           legend.position = "none",
#           strip.background = element_blank(),
#           strip.text = element_text(color = "grey10"),
#           text = element_text(size = 6)))
# 
# library(patchwork)
# 
# p1 + p2 + p3
# 
# ggsave("map_pagrus.png", width = 1280*2, height = 720*2, unit = "px")