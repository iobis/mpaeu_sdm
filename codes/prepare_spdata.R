library(obissdm)

# To track how it's going
st <- storr::storr_rds("datastd")

base_rast <- terra::rast("data/env/current/thetao_baseline_depthsurf_mean.tif")

obis <- read.csv("data/obis_splist_20230720.csv")
gbif <- read.csv("data/gbif_splist_20230720.csv")

library(dplyr)

obis <- obis %>% filter(kingdom != "Viruses" & kingdom != "Protozoa")
gbif <- gbif %>% filter(kingdom != "Viruses" & kingdom != "Protozoa")

sp_list <- c(obis$taxonID, gbif$valid_AphiaID)
sp_list <- unique(sp_list)

for (id in sp_list) {
  
  if (!st$exists(id)) {
    all_f <- list.files(paste0("data/species/key=", id), recursive = T)
    
    if (any(grepl("obis|gbif", all_f))) {
      mp_standardize(id,
                     geo_out_mccore = 4,
                     sdm_base = base_rast)
      st$set(id, "standardization_done")
    } else {
      st$set(id, "standardization_passed")
    }
    
  }
  cli::cli_progress_message("Concluded species {cli::bg_cyan(id)} ({cli::col_cyan(length(sp_list) - which(sp_list == id))} remaining)")
}

mp_get_ecoinfo(sp_list)
# 
# li <- ls()
# li <- li[!grepl("control|control_done", li)]
# rm(list = c(li, "li"))




for (id in sp_list) {
  
  if (!st$exists(id)) {
    
    all_f <- list.files(paste0("data/species/key=", id), recursive = T)
    
    if (any(grepl("dupclean", all_f))) {
      st$set(id, "standardization_done")
    }
  }
  
}


# If everything worked, remove storr
st$destroy()


for (i in 1:10) {
  cli::cli_progress_message("Tentando {i}")
  Sys.sleep(1)
}
