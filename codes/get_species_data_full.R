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
api_call <- httr::content(httr::GET("https://api.obis.org/export?complete=true"), as = "parsed")
api_call <- api_call$results[[1]]
latest_export <- paste0("https://obis-datasets.ams3.digitaloceanspaces.com/", api_call$s3path)

options(timeout=999999999)
download.file(url = latest_export,
              destfile = paste0("data/", gsub("exports/", "", api_call$s3path)),
              method = "wget")


# Get data from GBIF ----
sp_list_file <- list.files("data", pattern = "all_splist_", full.names = T)
sp_list <- read.csv(sp_list_file[order(sp_list_file, decreasing = T)][1])

# Because we already have the GBIF taxonKey we can supply this to the function
# Download and save
gbif_res <- mp_get_gbif(sci_names = na.omit(sp_list$gbif_speciesKey))

# It is likely that the download will fail, as our request is very big
# It can take up to 3h to process big requests
if (inherits(gbif_res, "occ_download")) {
    done <- FALSE
    trials <- 0
    while(!done & trials < 50) {
        if (interactive()) {
            if (!exists("timetry")) timetry <- as.numeric(readline("Chose time to try again, in minutes\n")) * 60
        } else {
            timetry <- 60 * 60
        }
        cat("Waiting... \n")
        Sys.sleep(timetry)
        cat("Trying to download\n")
        gbif_status <- try(rgbif::occ_download_meta(gbif_res)$status)
        if (gbif_status == "SUCCEEDED") {
            stop("deu")
            done <- TRUE
        } else if (gbif_status == "FAILED") {
            stop("Problem in the download.\n")
        } else {
            trials <- trials + 1
        }
    }
}

### END