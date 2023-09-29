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
library(biooracler)
library(cli)
library(terra)

# List datasets to download ----
datasets <- c(
  "thetao_baseline_2000_2019_depthsurf",
  "so_baseline_2000_2019_depthsurf",
  "PAR_mean_baseline_2000_2020_depthsurf",
  "po4_baseline_2000_2018_depthsurf",
  "phyc_baseline_2000_2020_depthsurf",
  "sws_baseline_2000_2019_depthsurf",
  "siconc_baseline_2000_2020_depthsurf"
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

# Define variables to be download
# In general, available are: min, mean, max, range, ltmin, ltmax, and sd
variables <- c("min", "mean", "max")

# Define out directory
fs::dir_create("data/raw/env_layers")
outdir <- "data/env/"
fs::dir_create(paste0(outdir, c("current", paste0("future/", future_scenarios), "terrain")))

# We create a function to check if the file is already downloaded
# If you want to ignore and download again, just set ignore = T
# in the loop functions
check_file_exist <- function(file_name, expres, ignore = FALSE) {
  if (!file.exists(file_name) | ignore) {
    expres
  } else {
    cli_alert_info("{file_name} already exists.")
  }
}


# Download all files ----

# We create a loop that will go through all the datasets
for (id in datasets) {
  cli_alert_info("Downloading {id}")
  
  # Generate the correct name for the variable, e.g. thethao_mean
  ds_vars <- gsub("_baseline(.*)$", "", id)
  ds_vars <- paste(ds_vars, variables, sep = "_")
  
  # Run for each time step
  for (ts in 1:length(time_steps)) {
    
    # Get name of the period (e.g. current)
    period <- names(time_steps)[ts]
    
    if (period == "current") {
      
      # Run for each variable
      for (sel_var in ds_vars) {
        
        # Get the final file name
        outfile <- paste0(outdir, "current/",
                          gsub("_2000_20[[:digit:]][[:digit:]]", "", tolower(id)),
                          "_", gsub("^(.*)_", "", sel_var), ".tif")
        
        # Check if exists
        check_file_exist(outfile, {
          
          # If it does not exists, try to download                        
          cli_progress_step("Downloading {sel_var} - {period}", spinner = T,
                            msg_failed = "Variable {id} [{sel_var}] not available")
          
          var <- try(download_dataset(tolower(id), sel_var, list(time = time_steps[[ts]]), fmt = "raster",
                                      directory = "data/raw/env_layers",
                                      verbose = F), silent = T) # Set verbose=TRUE to debug
          
          # If download succeed, save
          if (!assertthat::is.error(var)) {
            writeRaster(var, outfile, overwrite = T);rm(var)
          } else {
            cli_progress_done(result = "failed")
          }
          # To ignore the file checking and download anyway, set this to TRUE
       }, ignore = FALSE)
        cli_progress_done()
      }
      
    } else {
      
      # For the future we run for each scenario
      for (scen in future_scenarios) {
        # Get the modified dataset ID (each future have one ID)
        mod_id <- gsub("baseline_2000_20[[:digit:]][[:digit:]]",
                       paste0(scen, "_2020_2100"), id)
        
        # Run for each variable
        for (sel_var in ds_vars) {
          
          # Get the final file name
          outfile <- paste0(outdir, "future/", scen, "/",
                            gsub("_2020_2100", "", tolower(mod_id)),
                            "_", period, "_", gsub("^(.*)_", "", sel_var), ".tif")
          
          # Check if file exist
          check_file_exist(outfile, {
            
            # If it does not exists, try to download       
            cli_progress_step("Downloading {scen} - {sel_var} - {period}", spinner = T,
                              msg_failed = "Variable {mod_id}, scenario {scen} [{sel_var}], period {period} not available")
            
            var <- try(download_dataset(tolower(mod_id), sel_var, list(time = time_steps[[ts]]), fmt = "raster",
                                        directory = "data/raw/env_layers",
                                        verbose = F), silent = T) # Set verbose=TRUE to debug
            
            # If download succeed, save
            if (!assertthat::is.error(var)) {
              writeRaster(var, outfile, overwrite = T);rm(var)
            } else {
              cli_progress_done(result = "failed")
            }
            # To ignore the file checking and download anyway, set this to TRUE
          }, ignore = FALSE)
          cli_progress_done()
        }
      }
    }
  }
  # If everything is done, conclude last message
  cli_progress_done()
}


# Download terrain layers
terrain_vars <- c("bathymetry_mean")

for (tv in terrain_vars) {
  
  outfile <- paste0(outdir, "terrain/", tv, ".tif")
  
  check_file_exist(
    outfile,
    {
      # If it does not exists, try to download       
      cli_progress_step("Downloading terrain {tv}", spinner = T,
                        msg_failed = "Variable {tv} not available")
      
      var <- try(download_dataset("terrain_characteristics", tv,
                                  list(time = c("1970-01-01T00:00:00Z", "1970-01-01T00:00:00Z")), fmt = "raster",
                                  directory = "data/raw/env_layers",
                                  verbose = T), silent = T) # Set verbose=TRUE to debug
      
      # If download succeed, save
      if (!assertthat::is.error(var)) {
        writeRaster(var, outfile, overwrite = T);rm(var)
      } else {
        cli_progress_done(result = "failed")
      }
    }
  )
};cli_progress_done()


# Delete raw files (optional, but recommended) ----
fs::dir_delete("data/raw/env_layers")
