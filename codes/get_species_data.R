############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# June of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################### Download species data ##############################


# Load packages ----
library(robis)
library(rgbif)
library(obissdm)
library(arrow)
library(dplyr)


# Get data from OBIS ----
# Load species list
obis_splist <- read.csv("data/obis_splist_20230720.csv")

# Download and save
mp_get_obis(sci_names = obis_splist$taxonID, mode = "full")

# NOTE: for this specific project we used the full export of OBIS,
# which is available from https://obis.org/data/access/, thus, the data download
# step was skipped.


# Get data from GBIF ----
# Load species list
gbif_splist <- read.csv("data/gbif_splist_20231023.csv")

# Because we already have the GBIF taxonKey we can supply this to the function
# Download and save
mp_get_gbif(sci_names = gbif_splist$gbif_taxonKey)


# Split datasets ----
# Because we are working with a large quantity of species, we opt to split the
# dataset by key (taxon), so we can easily access each species without having to
# filter the full dataset each time

# Select columns to keep on files
gbif_columns <- c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                  "decimallongitude", "depth", "depthaccuracy", "day", "month",
                  "year", "taxonkey", "specieskey", "basisofrecord", "catalognumber")

obis_columns <- c("occurrenceID", "catalogNumber", "recordNumber", "fieldNumber", "AphiaID",
                  "materialSampleID", "institutionID", "collectionID", "datasetID",
                  "collectionCode", "institutionCode", "datasetName",
                  "eventID", "parentEventID", "decimalLatitude", "decimalLongitude", 
                  "species", "eventDate", "date_year", "day", "month", "year",
                  "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", 
                  "minimumDepthInMeters")

# Split GBIF dataset

# GBIF species are identified by a taxonKey, which differ from the AphiaID from 
# WoRMS used by OBIS. The function split_dataset enable to change the key for the
# new AphiaID key
gbif_new_keys <- gbif_splist[,c("gbif_speciesKey", "AphiaID")]
colnames(gbif_new_keys) <- c("key", "new_key")

split_dataset("data/raw/gbif_full_20230728.parquet",
              database_name = "gbif",
              grouping_key = "specieskey",
              sel_keys = gbif_new_keys$new_key,
              change_key = gbif_new_keys,
              sel_columns = gbif_columns,
              run_in_batches = T,
              batch_size = 100)

# Split OBIS dataset
split_dataset("data/raw/obis_20230726.parquet",
              database_name = "obis",
              grouping_key = "AphiaID",
              sel_keys = obis_splist$taxonID,
              sel_columns = obis_columns,
              run_in_batches = T,
              batch_size = 100)

# Check
sum(grepl("ftype=gbif", list.files("data/species", recursive = T)))
sum(grepl("ftype=obis", list.files("data/species", recursive = T)))

### END