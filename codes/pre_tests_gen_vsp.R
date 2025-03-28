############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############ SDM - generate virtual species for methods testing ################

# Load packages ----
library(obissdm)
library(terra)
source("functions/load_env.R")
set.seed(2023)



# List and load environmental layers ----
env_vars <- c("thetao-mean", "sws-mean", "so-mean", "o2-mean")

curr <- load_env(env_vars, terrain_vars = c("bathymetry_mean"))
ssp1 <- load_env(env_vars, scenario = "ssp126", terrain_vars = c("bathymetry_mean"))
ssp5 <- load_env(env_vars, scenario = "ssp585", terrain_vars = c("bathymetry_mean"))

# Load study area shapefile
starea <- vect("data/shapefiles/mpa_europe_starea_v3.gpkg")
exp_starea <- vect("data/shapefiles/mpa_europe_extarea_allatl_v3.gpkg")

# # Crop to the expanded area
# sampling_area <- crop(curr[[1]], exp_starea)
# sampling_area <- mask(sampling_area, exp_starea, inverse = F)
# sampling_area[!is.na(sampling_area)] <- 1
# sampling_area <- as.polygons(sampling_area)

# Create a bias layer that will give higher probability of sampling in shallow areas
# bias_layer <- curr$bathymetry_mean
# bias_layer[bias_layer >= -200] <- 1
# bias_layer[bias_layer < -200] <- 0.5
bias_layer <- rast("data/virtual_species/sampling_effort.tif")

# Load Gaussian field layer
gaus_rast <- rast("data/virtual_species/gaussian_field_layer.tif")


# Create Virtual Species ----

layers_1 <- curr
layers_2 <- c(curr, gaus_rast)

future_list <- list(
  ssp1 = ssp1, ssp5 = ssp5
)

future_list_2 <- list(
  ssp1 = c(ssp1, gaus_rast), ssp5 = c(ssp5, gaus_rast)
)

# VSP 1: coastal species, all variables known
vsp1_funs <- virtualspecies::formatFunctions(
  layers_1,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 6),
  sws_mean = c(fun = "dnorm", mean = 0.02, sd = 0.1),
  so_mean = c(fun = "dnorm", mean = 38, sd = 8),
  o2_mean = c(fun = "dnorm", mean = 400, sd = 100),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 150)
)

# VSP2: coastal species, not all variables known
vsp2_funs <- virtualspecies::formatFunctions(
  layers_2,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 6),
  sws_mean = c(fun = "dnorm", mean = 0.02, sd = 0.1),
  so_mean = c(fun = "dnorm", mean = 38, sd = 8),
  o2_mean = c(fun = "dnorm", mean = 400, sd = 100),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 150),
  random_gauss = c(fun = "dnorm", mean = 0.5, sd = 0.2)
)

gen_vsp(layers_1, vsp1_funs, vsp_name = "Virtual species 1",
        vsp_class = "all_var_coastal_vsp", save_key = "1001", plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list)

gen_vsp(layers_2, vsp2_funs, vsp_name = "Virtual species 2",
        vsp_class = "part_var_coastal_vsp", save_key = "1002", plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list_2)




# VSP 3: broad range species, all variables known
vsp3_funs <- virtualspecies::formatFunctions(
  layers_1,
  thetao_mean = c(fun = "dnorm", mean = 20, sd = 10),
  sws_mean = c(fun = "dnorm", mean = 0.02, sd = 0.1),
  so_mean = c(fun = "dnorm", mean = 30, sd = 10),
  o2_mean = c(fun = "dnorm", mean = 400, sd = 100),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 2000)
)

# VSP4: broad range species, not all variables known
vsp4_funs <- virtualspecies::formatFunctions(
  layers_2,
  thetao_mean = c(fun = "dnorm", mean = 20, sd = 10),
  sws_mean = c(fun = "dnorm", mean = 0.02, sd = 0.1),
  so_mean = c(fun = "dnorm", mean = 30, sd = 10),
  o2_mean = c(fun = "dnorm", mean = 400, sd = 100),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 2000),
  random_gauss = c(fun = "dnorm", mean = 0.5, sd = 0.2)
)

gen_vsp(layers_1, vsp3_funs, vsp_name = "Virtual species 3",
        vsp_class = "all_var_broad_vsp", save_key = "1003", plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = exp_starea,
        pred_fut = future_list)

gen_vsp(layers_2, vsp4_funs, vsp_name = "Virtual species 4",
        vsp_class = "part_var_broad_vsp", save_key = "1004", plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = exp_starea,
        pred_fut = future_list_2)


#END