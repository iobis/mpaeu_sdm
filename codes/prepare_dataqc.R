############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ Prepare occurrence data ###########################

# Previous step: see `get_species_data.R`

# Open packages ----
library(duckdb)
library(arrow)
library(dplyr)
library(storr)
library(furrr)


# Settings ----
shared <- "../mpaeu_shared/"

obis_dataset <- paste0(shared, "obis_20231025.parquet")
gbif_dataset <- paste0(shared, "gbif_20231025.parquet")

# Create database
db_file <- "data/species_qc.db"
con <- dbConnect(RSQLite::SQLite(), dbname = db_file)

# Create storr
st <- storr_rds("qc_storr")


# Load data ----
obis_pq <- open_dataset(obis_dataset) %>%
  select(all_of()) %>%
  filter() %>%
  to_duckdb()

gbif_pq <- open_dataset(gbif_dataset) %>%
  to_duckdb()

obis_list <- read.csv("data/obis_splist_20230720.csv")
gbif_list <- read.csv("data/gbif_splist_20231023.csv")

full_list <- bind_rows(gbif_list,
                       obis_list %>% rename(AphiaID = taxonID))

full_list <- full_list %>% distinct(AphiaID, .keep_all = T)


# Run QC in parallel ----
plan(multisession, workers = 5)

future_map(1:nrow(full_list), function(id){
  
  obis_id <- full_list$AphiaID[id]
  gbif_id <- full_list$gbif_taxonKey[id]
  
  obis_data <- obis_pq %>%
    filter(AphiaID == obis_id) %>%
    collect()
  
  if (!is.na(gbif_id)) {
    gbif_data <- gbif_pq %>%
      filter(AphiaID == gbif_id) %>%
      collect() # TODO: Rename any column needed!!!!!!
  } else {
    gbif_data <- NULL
  }
  
  full_data <- bind_rows(obis_data, gbif_data)
  
  std_data <- obissdm::mp_standardize(full_data)
  
  if (!is.null(std_data)) {
    st$set(id, std_data)
    to_return <- "done"
  } else {
    to_return <- "failed"
  }
  
  return(to_return)
  
}, .progress = T)