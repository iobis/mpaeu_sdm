############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
###################### Download environmental layers ###########################

# Environmental layers are downloaded from the Bio-ORACLE ERDDAP server
# To see all available options, visit:
# https://erddap.bio-oracle.org/erddap/griddap/index.html?page=1&itemsPerPage=1000
# We use the biooracler package which is still not on CRAN (but soon will be)
# The files are saved on the following structure
# data
# \_ raw
#.     \_ env_layers: raw files downloaded from Bio-ORACLE - deleted at the end
# \_ env: environmental layers
#      \_ terrain: terrain layers (static)
#.     \_ current: current period layers
#      \_ future: future period scenarios folders
#.          \_ ssp1: future layers for SSP1 scenario
#           \_ ...: other scenarios


# Load packages ----
#devtools::install_github("bio-oracle/biooracler")
library(obissdm)


# List datasets to download ----
datasets <- c(
  "thetao_baseline_2000_2019_depthsurf",
  "so_baseline_2000_2019_depthsurf",
  "PAR_mean_baseline_2000_2020_depthsurf",
  "phyc_baseline_2000_2020_depthsurf",
  "ph_baseline_2000_2018_depthsurf",
  "sws_baseline_2000_2019_depthsurf",
  "siconc_baseline_2000_2020_depthsurf",
  "o2_baseline_2000_2018_depthsurf",
  "KDPAR_mean_baseline_2000_2020_depthsurf",
  "no3_baseline_2000_2018_depthsurf",
  "chl_baseline_2000_2018_depthsurf",
  "tas_baseline_2000_2020_depthsurf",
  "si_baseline_2000_2018_depthsurf",
  "mlotst_baseline_2000_2019_depthsurf"
)

datasets <- c(datasets,
              gsub("depthsurf", "depthmean", datasets))


# List scenarios to download ----
future_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585")


# Define time steps ----
time_steps <- list(
  current = c("2000-01-01T00:00:00Z", "2010-01-01T00:00:00Z"), #2000-2010/2010-2020
  dec50 = c("2030-01-01", "2040-01-01"), #2030-2040/2040-2050
  dec100 = c("2080-01-01", "2090-01-01") #2080-2090/2090-2100
)

# Define variables to be downloaded
# In general, available are: min, mean, max, range, ltmin, ltmax, and sd
variables <- c("min", "mean", "max")


get_env_data(datasets = datasets, future_scenarios = future_scenarios,
             time_steps = time_steps, variables = variables,
             terrain_vars = c(
               "bathymetry_mean",
               "slope",
               "terrain_ruggedness_index"
             ), average_time = T)

# For just temperature, download also the range, ltmin and ltmax
get_env_data(datasets = "thetao_baseline_2000_2019_depthsurf",
             future_scenarios = future_scenarios,
             time_steps = time_steps, variables = c("range", "ltmin", "ltmax"),
             average_time = T)

# For Chlorophyll-a we remove the depthmean and depthmax, as for the future is
# not available
to_remove <- list.files("data/env/current", full.names = T)
to_remove <- to_remove[grepl("chl", to_remove)]
to_remove <- to_remove[grepl("depthmean|depthmax", to_remove)]
fs::file_delete(to_remove)

# Rename KDPAR for kd
to_rename <- list.files("data/env", recursive = T, full.names = T)
to_rename <- to_rename[grepl("kdpar", to_rename)]
new_names <- gsub("kdpar", "kd", to_rename)
file.rename(to_rename, new_names)

# Rename terrain_ruggedness
to_rename <- list.files("data/env/terrain/", recursive = T, full.names = T)
to_rename <- to_rename[grepl("rugg", to_rename)]
to_rename <- to_rename[!grepl("aux", to_rename)]
new_names <- gsub("terrain_ruggedness_index", "rugosity", to_rename)
edit_r <- terra::rast(to_rename)
names(edit_r) <- "rugosity"
terra::writeRaster(edit_r, new_names, overwrite = T)
fs::file_delete(to_rename)
### END