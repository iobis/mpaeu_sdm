############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################## Create additional environmental layers ######################

# Load packages ----
library(terra)
library(stars)


# Wavefetch aggregated layers ----
download.file("https://figshare.com/ndownloader/articles/8668127/versions/1",
  destfile = "data/raw/wavefetch.zip",
  method = "wget")

unzip("data/raw/wavefetch.zip", exdir = "data/raw/wavefetch")
file.remove("data/raw/wavefetch.zip")

lf <- list.files("data/raw/wavefetch", full.names = T)
lf <- lf[!grepl("aux", lf)]

mosaic <- eval(parse(text = paste(
  "st_mosaic(", paste0("read_stars('", lf , "')", collapse = ","), ")"
)))

plot(mosaic)

write_stars(mosaic, "data/temp_wavef.tif")

mosaic_rast <- rast("data/temp_wavef.tif")

base_file <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")

europe <- vect("data/shapefiles/mpa_europe_starea_v2.shp")

reproj <- project(mosaic_rast, base_file)
reproj <- crop(reproj, ext(europe))

plot(reproj)
lines(europe)

names(reproj) <- "wavefetch"
writeRaster(reproj, "data/env/terrain/wavefetch.tif", overwrite = T)

fs::file_delete("data/temp_wavef.tif")


# Distance to coast layer ----
# Load a base layer
base <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")

coast <- base
coast[] <- NA

coast_mask <- mask(coast, base, updatevalue = 1, inverse = F)
plot(coast_mask)

coast_agg <- aggregate(coast_mask, 4, na.rm = T)

coast_dist <- terra::distance(coast_agg)
plot(coast_dist)

coast_dist <- disagg(coast_dist, 4)
coast_dist <- mask(coast_dist, base)

coast_dist <- coast_dist/1000 # to km

names(coast_dist) <- "coastdist"

writeRaster(coast_dist, "data/env/terrain/distcoast.tif", overwrite = T)



# Distance to closest estuary ----
# Load estuaries
# Data for estuaries was downloaded from here: https://data.unep-wcmc.org/datasets/23
estuary_points <- vect("data/raw/estuaries2003_v2/14_001_UBC003_SAU_Estuaries2003_v2.shp")

base <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")

est_rast <- base
est_rast[!is.na(est_rast)] <- 0

est_rast_agg <- aggregate(est_rast, 2, na.rm = T)

est_rast_agg[estuary_points] <- 1

distestuary_mean <- terra::gridDist(est_rast_agg, target = 1, scale = 1000)

distestuary_mean <- disagg(distestuary_mean, 2)
distestuary_mean <- mask(distestuary_mean, base)
plot(distestuary_mean)

names(distestuary_mean) <- "distestuary"

writeRaster(distestuary_mean, "data/env/terrain/distestuary.tif", overwrite = T)



# Distance to sea icea ----
all_siconc <- list.files("data/env", recursive = T, full.names = T)
all_siconc <- all_siconc[grepl("siconc", all_siconc)]
all_siconc <- all_siconc[!grepl("aux", all_siconc)]
all_siconc <- all_siconc[grepl("_mean", all_siconc)]

for (to_do in all_siconc) {
  
  cat("Processing", to_do, "\n")
  
  # Load sea ice layer
  seaice <- rast(to_do)
  
  seaice_agg <- aggregate(seaice, 2, na.rm = T)
  
  seaice_agg[seaice_agg > 0] <- 1 
  
  seaicedist_mean <- terra::gridDist(seaice_agg, target = 1, scale = 1000)
  
  seaicedist_mean <- disagg(seaicedist_mean, 2)
  
  seaicedist_mean <- mask(seaicedist_mean, seaice)
  
  #plot(seaicedist_mean)
  
  names(seaicedist_mean) <- "seaicedist"
  
  outfile <- gsub("siconc", "seaicedist", to_do)
  
  writeRaster(seaicedist_mean, outfile, overwrite = T)
}
