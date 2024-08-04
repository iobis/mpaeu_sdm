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
destroy_st <- FALSE

# Create storr
st <- storr_rds("qc_storr")

# Open target species list
sp_list_p <- recent_file("data", "all_splist")
species_list <- read.csv(sp_list_p)

priority <- read.csv("../mpaeu_constatus/v4/MPAEU_D3_4_threatened_species_fulltaxonomy.csv")
priority <- priority$valid_AphiaID
priority <- priority[!is.na(priority)]
species_list <- species_list[species_list$AphiaID %in% priority,]

# Run QC in parallel ----
plan(multisession, workers = 3) # Need to be multisession, otherwise crashes

proc_res <- future_map(1:nrow(species_list), function(id){
  
  require(terra)

  sp <- species_list$taxonID[id]
  
  if (!st$exists(sp)) {
    
    res <- try(mp_standardize(species = sp,
                              sdm_base = terra::rast("data/env/current/thetao_baseline_depthsurf_mean.tif"),
                              species_list = sp_list_p,
                              species_folder = "data/raw/",
                              reader = "arrow"))
    
    if (!inherits(res, "try-error")) {
      st$set(sp, "done")
      to_return <- "done"
    } else {
      to_return <- list("failed", res)
      st$set(sp, to_return)
    }
  } else {
    to_return <- "already_done"
  }
  
  return(to_return)
  
}, .progress = T)

# Check if everything was done
if (length(species_list$taxonID[!species_list$taxonID %in% as.numeric(st$list())]) > 0) {
  stop("Not all species were processed. Check.")
} else {
  cli::cli_alert_success("All species processed.")
  if (destroy_st) {
    # Destroy st
    st$destroy()
  }
}

### END