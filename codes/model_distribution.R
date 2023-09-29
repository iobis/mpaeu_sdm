############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### Model the distribution of species ######################


# Load packages / settings ----
library(obissdm)
library(arrow)
library(terra)
library(parallel)
source("functions/load_env.R")


# Load species list ----


# Load list of environmental variables ----
env_vars <- c("thetao-mean", "so-mean", "po4-mean", "phyc-mean")

curr <- load_env(env_vars, terrain_vars = "bathymetry_mean")
ssp1 <- load_env(env_vars, scenario = "ssp126", terrain_vars = "bathymetry_mean")
ssp2 <- load_env(env_vars, scenario = "ssp245", terrain_vars = "bathymetry_mean")
ssp3 <- load_env(env_vars, scenario = "ssp370", terrain_vars = "bathymetry_mean")
ssp4 <- load_env(env_vars, scenario = "ssp460", terrain_vars = "bathymetry_mean")
ssp5 <- load_env(env_vars, scenario = "ssp585", terrain_vars = "bathymetry_mean")

# Load study area shapefile
starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
exp_starea <-  ext(-41, 47, 20, 89) # ext(starea) +- 5 degrees

# Crop to the expanded area (only the one that will be used for
# sampling the background)

# Remove Red sea
mregions <- mregions::mr_shp("MarineRegions:iho")
mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

curr <- mask(crop(curr, exp_starea), mregions, inverse = T)

# Mask to a depth of up to 1500 m
bath_mask <- curr$bathymetry_mean
bath_mask[bath_mask < -1500] <- NA

curr <- mask(curr, bath_mask)

ssp1 <- mask(crop(ssp1, exp_starea), curr)
ssp2 <- mask(crop(ssp2, exp_starea), curr)
ssp3 <- mask(crop(ssp3, exp_starea), curr)
ssp4 <- mask(crop(ssp4, exp_starea), curr)
ssp5 <- mask(crop(ssp5, exp_starea), curr)

block_list <- list(
  spatial_grid = rast(ext(curr), resolution = 5),
  spatial_lat = rast(ext(curr), ncol = 1, nrow = 30)
)


# Load list of mode of life ----
species_list <- c(127053, 127054, 127063)


# Create the function with all the steps ----
model_sp <- function(species
                     #, env_list, life_mode
                     ) {
  
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
                             env_layers = scale(curr),
                             quad_number = 150000)
  
  sp_data <- mp_prepare_blocks(sp_data, method = "manual", 
                               manual_shp = block_list,
                               n_iterate = 300)
  
  # Model
  sp_model <- sdm_module_lasso(sp_data, blocks_all = TRUE)
  
  # Predict
  c_m <- global(curr, mean, na.rm = T)[,1]
  c_sd <- global(curr, sd, na.rm = T)[,1]
  curr <- (curr - c_m) / c_sd
  ssp1 <- (ssp1 - c_m) / c_sd
  ssp2 <- (ssp2 - c_m) / c_sd
  ssp3 <- (ssp3 - c_m) / c_sd
  ssp4 <- (ssp4 - c_m) / c_sd
  ssp5 <- (ssp5 - c_m) / c_sd
  
  curr_pred <- predict(sp_model, curr)
  ssp1_pred <- predict(sp_model, ssp1)
  ssp2_pred <- predict(sp_model, ssp2)
  ssp3_pred <- predict(sp_model, ssp3)
  ssp4_pred <- predict(sp_model, ssp4)
  ssp5_pred <- predict(sp_model, ssp5)
  
  mod_resp_cur <- resp_curves(sp_model, scale(curr))
  
  return(list(
    current = curr_pred,
    ssp1 = ssp1_pred,
    ssp2 = ssp2_pred,
    ssp3 = ssp3_pred,
    ssp4 = ssp4_pred,
    ssp5 = ssp5_pred,
    resp_cur = mod_resp_cur,
    pts = pts
  ))
  
}


# Run in parallel ----
pre_results <- mclapply(species_list, model_sp, mc.cores = 3)

teste <- model_sp(species_list[1])
