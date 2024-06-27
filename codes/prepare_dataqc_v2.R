############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ Prepare occurrence data ###########################

# Previous step: see `get_species_data.R`

# Load packages ----
library(obissdm)
library(terra)
library(storr)
library(furrr)

# Create storr
st <- storr_rds("qc_storr")

# Open target species list
species_list <- read.csv("analysis/internal_eval_2024/species_list.csv")

# Run QC in parallel ----
plan(multicore, workers = 6)

proc_res <- future_map(1:nrow(species_list), function(id){
  
  #library(terra)
  
  sp <- species_list$taxonID[id]
  
  if (!st$exists(sp)) {
    
    res <- try(mp_standardize(species = sp,
                              sdm_base = rast("data/env/current/thetao_baseline_depthsurf_mean.tif"),
                              species_list = "data/all_splist_20240319.csv",
                              species_folder = "../mpaeu_shared/"))
    
    if (!inherits(res, "try-error")) {
      st$set(sp, "done")
      to_return <- "done"
    } else {
      to_return <- list("failed", res)
    }
  } else {
    to_return <- "already_done"
  }
  
  return(to_return)
  
}, .progress = T)

# Check if everything was done
species_list$taxonID[!species_list$taxonID %in% st$list()]

# Destroy st
st$destroy()
