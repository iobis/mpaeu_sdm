library(obissdm)
library(arrow)

gbif_list <- read.csv("data/gbif_splist_20230720.csv")
obis_list <- read.csv("data/obis_splist_20230720.csv")

mp_get_local(
  "data/raw/gbif_full_20230728.parquet",
  database_name = "gbif",
  gbif_key = gbif_list$gbif_taxonKey,
  aphia_id = gbif_list$AphiaID,
  sel_columns = c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                  "decimallongitude", "depth", "depthaccuracy", "day", "month",
                  "year", "taxonkey", "basisofrecord", "catalognumber")
)

get_conf <- readline("Should continue? y/n")

if (get_conf == "y") {
  mp_get_local(
    "data/raw/obis_20230726.parquet",
    database_name = "obis",
    aphia_id = obis_list$taxonID,
    db_mode = TRUE,
    sel_columns = c("occurrenceID", "catalogNumber", "recordNumber", "fieldNumber",
                    "materialSampleID", "institutionID", "collectionID", "datasetID",
                    "collectionCode", "institutionCode", "datasetName", "eventID", "parentEventID",
                    "decimalLatitude", "decimalLongitude", "species", "eventDate", "date_year", "day", "month", "year",
                    "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", "minimumDepthInMeters")
  )
}


bench::mark(
  {mp_get_local(
    "data/raw/gbif_full_20230728.parquet",
    database_name = "gbif",
    gbif_key = gbif_list$gbif_taxonKey[1:5],
    aphia_id = gbif_list$AphiaID[1:5],
    sel_columns = c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                    "decimallongitude", "depth", "depthaccuracy", "day", "month",
                    "year", "taxonkey", "basisofrecord", "catalognumber")
  )},
  {mp_get_local(
    "data/raw/gbif_full_20230728.parquet",
    database_name = "gbif",
    gbif_key = gbif_list$gbif_taxonKey[1:5],
    aphia_id = gbif_list$AphiaID[1:5],
    sel_columns = c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                    "decimallongitude", "depth", "depthaccuracy", "day", "month",
                    "year", "taxonkey", "basisofrecord", "catalognumber"),
    db_mode = T
  )},
  check = F
)
