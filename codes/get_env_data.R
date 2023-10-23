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
  "po4_baseline_2000_2018_depthsurf",
  "phyc_baseline_2000_2020_depthsurf",
  "ph_baseline_2000_2018_depthsurf",
  "sws_baseline_2000_2019_depthsurf",
  "siconc_baseline_2000_2020_depthsurf",
  "o2_baseline_2000_2018_depthsurf",
  "KDPAR_mean_baseline_2000_2020_depthsurf",
  "dfe_baseline_2000_2018_depthsurf",
  "no3_baseline_2000_2018_depthsurf",
  "chl_baseline_2000_2018_depthsurf",
  "tas_baseline_2000_2020_depthsurf"
)

datasets <- c(datasets,
              gsub("depthsurf", "depthmean", datasets),
              gsub("depthsurf", "depthmax", datasets))

# List scenarios to download ----
future_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585")

# Define time steps ----
time_steps <- list(
  current = c("2010-01-01", "2010-01-01"),
  dec50 = c("2050-01-01", "2050-01-01"),
  dec100 = c("2090-01-01", "2090-01-01")
)

# Define variables to be downloaded
# In general, available are: min, mean, max, range, ltmin, ltmax, and sd
variables <- c("min", "mean", "max")


get_env_data(datasets = datasets, future_scenarios = future_scenarios,
             time_steps = time_steps, variables = variables,
             terrain_vars = "bathymetry_mean")

# For just temperature, download also the range
get_env_data(datasets = "thetao_baseline_2000_2019_depthsurf",
             future_scenarios = future_scenarios,
             time_steps = time_steps, variables = c("range"))
