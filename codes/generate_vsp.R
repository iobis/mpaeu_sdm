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
library(raster) # needed because virtualspecies package is still using raster
source("functions/load_env.R")
set.seed(2023)



# List and load environmental layers ----
env_vars <- c("thetao-mean", "so-mean", "po4-mean", "phyc-mean")

curr <- load_env(env_vars, terrain_vars = "bathymetry_mean")
ssp1 <- load_env(env_vars, scenario = "ssp126", terrain_vars = "bathymetry_mean")
ssp5 <- load_env(env_vars, scenario = "ssp585", terrain_vars = "bathymetry_mean")

# Load study area shapefile
starea <- shapefile("data/shapefiles/mpa_europe_starea_v2.shp")
exp_starea <- shapefile("data/shapefiles/mpa_europe_extarea_allatl_v2.shp")

# Crop to the expanded area
curr <- crop(curr, exp_starea)
ssp1 <- crop(ssp1, exp_starea)
ssp5 <- crop(ssp5, exp_starea)

# Remove Pacific and Indian ocean
mregions <- mregions::mr_shp("MarineRegions:iho")
mregions <- mregions[mregions$name %in% 
                       c("South Pacific Ocean", "North Pacific Ocean",
                         "Indian Ocean", "Mozambique Channel", "Red Sea"),]

curr <- mask(curr, mregions, inverse = T)
ssp1 <- mask(ssp1, mregions, inverse = T)
ssp5 <- mask(ssp5, mregions, inverse = T)

# Create a bias layer that will give higher probability of sampling in shallow areas
bias_layer <- curr$bathymetry_mean
bias_layer[bias_layer >= -200] <- 1
bias_layer[bias_layer < -200] <- 0.5
# bias_layer <- (bias_layer - global(bias_layer, min, na.rm = T)[,1]) /
#   (global(bias_layer, max, na.rm = T)[,1] - global(bias_layer, min, na.rm = T)[,1])
# bias_layer[bias_layer < 0.1] <- 0.1
bias_layer <- raster(bias_layer)

# Mask everything to a depth of up to 1500 m
bath_mask <- curr$bathymetry_mean
bath_mask[bath_mask < -1500] <- NA

curr <- mask(curr, bath_mask)
ssp1 <- mask(ssp1, bath_mask)
ssp5 <- mask(ssp5, bath_mask)



# Create Virtual Species ----

# VSP 1 and 2 (all variables known, broad distribution vs Europe limited)
# Convert to raster format
layers_1 <- stack(curr)
layers_ssp1 <- stack(ssp1)
layers_ssp5 <- stack(ssp5)
future_list <- list(ssp1 = layers_ssp1,
                    ssp5 = layers_ssp5)

# Get functions
vsp1_funs <- virtualspecies::formatFunctions(
  layers_1,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 6),
  so_mean = c(fun = "dnorm", mean = 38, sd = 8),
  po4_mean = c(fun = "dnorm", mean = 0, sd = 0.7),
  phyc_mean = c(fun = "dnorm", mean = 2.5, sd = 3.5),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 400)
)

vsp2_funs <- virtualspecies::formatFunctions(
  layers_1,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 2),
  so_mean = c(fun = "dnorm", mean = 38, sd = 6),
  po4_mean = c(fun = "dnorm", mean = 0, sd = 0.7),
  phyc_mean = c(fun = "dnorm", mean = 2.5, sd = 3.5),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 300)
)

gen_vsp(layers_1, vsp1_funs, vsp_name = "Virtual species 1",
        vsp_class = "all_var_wide_vsp", save_key = 101, plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list)


# VSP 2
# Same configuration as vsp 1, but constraining to Europe
gen_vsp(layers_1, vsp2_funs, vsp_name = "Virtual species 2",
        vsp_class = "all_var_narrow_vsp", save_key = 102, plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list)



# VSP 3 and 4 (unknown variable, broad distribution vs Europe limited)
# Generate a random Gaussian field that would represent a variable that
# is not being considered in the modeling process
random_field <- prioritizr::simulate_species(curr$thetao_mean, n = 1, scale = 0.1)
names(random_field) <- "gau"

layers_2 <- c(curr, random_field)

layers_2 <- stack(layers_2)

layers_ssp1 <- stack(c(ssp1, random_field))
layers_ssp5 <- stack(c(ssp5, random_field))
future_list <- list(ssp1 = layers_ssp1,
                    ssp5 = layers_ssp5)
linfun <- function(x, a, b) {(a*x)+b}

# Get functions
vsp3_funs <- virtualspecies::formatFunctions(
  layers_2,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 6),
  so_mean = c(fun = "dnorm", mean = 38, sd = 8),
  po4_mean = c(fun = "dnorm", mean = 0, sd = 0.7),
  phyc_mean = c(fun = "dnorm", mean = 2.5, sd = 3.5),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 400),
  gau = c(fun = "dnorm", mean = 1, sd = 1)
)

vsp4_funs <- virtualspecies::formatFunctions(
  layers_2,
  thetao_mean = c(fun = "dnorm", mean = 12, sd = 2),
  so_mean = c(fun = "dnorm", mean = 38, sd = 6),
  po4_mean = c(fun = "dnorm", mean = 0, sd = 0.7),
  phyc_mean = c(fun = "dnorm", mean = 2.5, sd = 3.5),
  bathymetry_mean = c(fun = "dnorm", mean = 0, sd = 300),
  gau = c(fun = "dnorm", mean = 1, sd = 1)
)

# VSP 3
gen_vsp(layers_2, vsp3_funs, vsp_name = "Virtual species 3",
        vsp_class = "gau_var_wide_vsp", save_key = 103, plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list)


# VSP 4
# Same configuration as vsp 1, but constraining to Europe
gen_vsp(layers_2, vsp4_funs, vsp_name = "Virtual species 4",
        vsp_class = "gau_var_narrow_vsp", save_key = 104, plot = F,
        samp_bias_layer = bias_layer, samp_constr_shape = starea,
        pred_fut = future_list)

# Save random field raster
writeRaster(layers_2$gau, "data/virtual_species/gaussian_field.tif")
