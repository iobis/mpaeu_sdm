library(obissdm)
library(arrow)
library(dplyr)

gbif_list <- read.csv("data/gbif_splist_20230720.csv")
obis_list <- read.csv("data/obis_splist_20230720.csv")

gbif_columns <- c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                  "decimallongitude", "depth", "depthaccuracy", "day", "month",
                  "year", "taxonkey", "specieskey", "basisofrecord", "catalognumber")

obis_columns <- c("occurrenceID", "catalogNumber", "recordNumber", "fieldNumber", "AphiaID",
                "materialSampleID", "institutionID", "collectionID", "datasetID",
                "collectionCode", "institutionCode", "datasetName", "eventID", "parentEventID",
                "decimalLatitude", "decimalLongitude", "species", "eventDate", "date_year", "day", "month", "year",
                "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", "minimumDepthInMeters")


### GBIF
gbif_new_keys <- gbif_list[,c("gbif_speciesKey", "AphiaID")]
colnames(gbif_new_keys) <- c("key", "new_key")

split_dataset("data/raw/gbif_full_20230728.parquet",
              database_name = "gbif",
              grouping_key = "specieskey",
              sel_keys = gbif_new_keys$new_key,
              change_key = gbif_new_keys,
              sel_columns = gbif_columns,
              run_in_batches = T,
              batch_size = 100)

# Rename decimalLatitude and decimalLongitude columns in GBIF datasets
all_f <- list.files("data/species/", recursive = T, full.names = T)
all_f <- all_f[grep("gbif", all_f)]

lapply(all_f, function(x){
  rd <- read_parquet(x)
  rd <- rd %>%
    rename(decimalLatitude = decimallatitude,
           decimalLongitude = decimallongitude)
  write_parquet(rd, x)
  return(invisible(NULL))
})

### OBIS
split_dataset("data/raw/obis_20230726.parquet",
              database_name = "obis",
              grouping_key = "AphiaID",
              sel_keys = obis_list$taxonID,
              sel_columns = obis_columns,
              run_in_batches = T,
              batch_size = 100)
