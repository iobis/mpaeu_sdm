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
fs::dir_create("data/raw/")

# Get data from OBIS ----
api_call <- httr::content(httr::GET("https://api.obis.org/export?complete=true"), as = "parsed")
api_call <- api_call$results[[1]]
latest_export <- paste0("https://obis-datasets.ams3.digitaloceanspaces.com/", api_call$s3path)
# NOTE: deprecated on January 2026 - use open data pathway now: https://github.com/iobis/obis-open-data

options(timeout = 999999999)
download.file(
    url = latest_export,
    destfile = paste0("data/raw/", gsub("exports/", "", api_call$s3path)),
    method = "wget"
)


# Get data from GBIF ----
sp_list_file <- recent_file("data", "all_splist")
sp_list <- read.csv(sp_list_file)

# In the case of our project, the data to be downloaded is very big. 
# There are 3 ways of downloading, each with benefits and drawbacks
# 1. download of subset through AWS: this will use the most recent monthly AWS 
# export of GBIF data. It uses the DuckDB backend to subset and download/save
# the data locally. It can take more than 3h, but return only the subset (we
# use this method here)
# 2. download full AWS export and subset: this may be faster if there is much data
# but needs more than 250GB of space available.
# 3. use the GBIF download API: takes more than 3h and also download a big file
# (because it includes media/issues). However, ensure that you have the most recent
# data (AWS exports are monthly)
# Codes for the three methods are below, but we use method 1.


# Download of subset through AWS ----
# Estimated time is >3h
gbif_res <- occurrence_gbif_db(full_mode = F, export_path = paste0("data/raw/gbif_", format(Sys.Date(), "%Y%m%d")),
                   scientificname = sp_list$gbif_speciesKey[!is.na(sp_list$gbif_speciesKey)],
                   verbose = T)

if(any(grepl("order=NULL", list.files(gbif_res)))) {
  file.rename(paste0(gbif_res, "/order=NULL"),
             paste0(gbif_res, "/order=noorder"))
}


# # Download full export and subset ----
# gbif_res <- occurrence_gbif_db(export_path = "data/raw", verbose = T)

# gbif_ds <- open_dataset(gbif_res)

# gbif_ds %>%
#   select(-identifiedby, -recordedby, -typestatus, -mediatype, -issue) %>%
#   filter(taxonkey %in% na.omit(sp_list$gbif_speciesKey)) %>%
#   group_by(order) %>%
#   write_dataset(path = paste0("gbif_", format(Sys.Date(), "%Y%m%d")),
#                 format = "parquet")

# file.rename(paste0("gbif_", format(Sys.Date(), "%Y%m%d"), "/order=__HIVE_DEFAULT_PARTITION__"),
#             paste0("gbif_", format(Sys.Date(), "%Y%m%d"), "/order=noorder"))

# # Optional: remove the full export
# # fs::file_delete(gbif_res)


# Download through the GBIF API ----
# # Because we already have the GBIF taxonKey we can supply this to the function
# # Download and save
# gbif_res <- mp_get_gbif(sci_names = na.omit(sp_list$gbif_speciesKey))

# # It is likely that the download will fail, as our request is very big
# # It can take up to 3h to process big requests
# if (inherits(gbif_res, "occ_download")) {
#     done <- FALSE
#     trials <- 0
#     while (!done & trials < 50) {
#         if (interactive()) {
#             if (!exists("timetry")) timetry <- as.numeric(readline("Chose time to try again, in minutes\n")) * 60
#         } else {
#             timetry <- 60 * 60
#         }
#         cat("Waiting... \n")
#         Sys.sleep(timetry)
#         cat("Trying to download\n")
#         gbif_status <- try(rgbif::occ_download_meta(gbif_res)$status)
#         if (gbif_status == "SUCCEEDED") {
#             dl <- rgbif::occ_download_get(gbif_res, path = "data/raw", curlopts = list(timeout_ms = 10000))
#             unzip(paste0("data/raw/", gbif_res, ".zip"),
#                 exdir = "data/raw/"
#             )

#             file.rename(
#                 "data/raw/occurrence.parquet",
#                 glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.parquet")
#             )
#             fs::file_delete(paste0(
#                 "data/raw/", gbif_res,
#                 ".zip"
#             ))

#             write.table(
#                 data.frame(
#                     date = attr(
#                         gbif_res,
#                         "created"
#                     ), download_code = as.character(gbif_res),
#                     doi = attr(gbif_res, "doi"), citation = attr(
#                         gbif_res,
#                         "citation"
#                     )
#                 ), glue::glue("data/gbif_full_download_log.txt"),
#                 row.names = F
#             )

#             done <- TRUE
#         } else if (gbif_status == "FAILED") {
#             stop("Problem in the download.\n")
#         } else {
#             trials <- trials + 1
#         }
#     }
# }

### END