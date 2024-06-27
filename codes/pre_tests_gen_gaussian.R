############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############ SDM - generate virtual species for methods testing ################

# Load packages ----
library(terra)
library(INLA)
library(inlabru)
library(fmesher)
library(sp)
set.seed(2023)

# Create Random Gaussian field layer
boundary <- spoly(data.frame(easting = c(-180, 180, 180, -180), northing = c(-90, -90, 90, 90)))
mesh_fine <- fm_mesh_2d_inla(boundary = boundary, max.edge = 3)

# Create matern field
matern_fine <-
  inla.spde2.pcmatern(mesh_fine,
                      prior.sigma = c(1, 0.01),
                      prior.range = c(1, 0.01)
  )
# Set parameters
true_range <- 80
true_sigma <- 5
true_q <- inla.spde.precision(matern_fine, theta = log(c(true_range, true_sigma)))

true_sd <- diag(inla.qinv(true_q)) ^ 0.5

# True field
true_field <- inla.qsample(1, true_q)[, 1]

truth <- expand.grid(
  easting = seq(-180, 180, length = 1800),
  northing = seq(-90, 90, length = 900)
)
truth <- sf::st_as_sf(truth, coords = c("easting", "northing"))
truth$field <- fm_evaluate(
  mesh_fine,
  loc = truth,
  field = true_field
)

# Convert to raster
gaus_rast <- rast(ncol = 1800, nrow = 900)
values(gaus_rast) <- truth$field

# Put in the same resolution as the layers
gaus_rast <- disagg(gaus_rast, fact = 4)

# Normalize to 0-1 scale
gaus_rast <- (gaus_rast - minmax(gaus_rast)[1,])/(minmax(gaus_rast)[2,] - minmax(gaus_rast)[1,])

# Mask according to the other layers
curr <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")
gaus_rast <- mask(gaus_rast, curr$thetao_mean)
names(gaus_rast) <- "random_gauss"
plot(gaus_rast)

# Save for later use
writeRaster(gaus_rast, "data/virtual_species/gaussian_field_layer.tif")
