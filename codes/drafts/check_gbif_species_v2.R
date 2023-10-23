#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe project
# s.principe@unesco.org

#### Get list of marine species occurring on the study area on GBIF

# Load packages
library(rgbif)
library(sf)
library(terra)
library(tidyverse)
library(worrms)
library(cli)
sf_use_s2(FALSE)

# Load study area
starea <- st_read("data/shapefiles/mpa_europe_starea_v2.shp")
starea <- st_buffer(starea, 0.5)

# We first download a list of all species occurring in a bounding box around the
# whole study area. We do this because our study area is complex, and is not
# possible to parse the geometry to the search.
splist_download <- occ_download(pred("occurrenceStatus", "PRESENT"),
                                pred_within(st_as_text(st_as_sfc(st_bbox(starea)))),
                                format = "SPECIES_LIST")

# Check if download is ready
occ_download_wait(splist_download)

# If ready, download
splist_down_get <- occ_download_get(splist_download)

# Import and delete zip file
splist_bbox <- occ_download_import(splist_down_get)
file.remove(splist_down_get)

# We then obtain the list of species only on the study area
# using the function occ_count - because the function will return only
# the taxon keys, we will need the list we obtained earlier to do the match
# (otherwise retrieving the names using name_usage for example would take long)

# Because the study area is complex, we divide it in blocks to run
# each one separately
grid <- st_make_grid(starea, n = c(25, 25))

starea_grid <- st_intersection(starea, grid)

# Get information from GBIF
# Set a progress bar to track results
cli_progress_bar(
  format = paste0(
    "{pb_spin} Getting species list - Number of species until now: {nrow(counts_list)} ",
    "{pb_bar} {pb_percent} ETA:{pb_eta}"
  ),
  format_done = paste0(
    "{col_green(symbol$tick)} Retrieved {nrow(counts_list)} species ",
    "in {pb_elapsed}."
  ),
  total = nrow(starea_grid)
)
for (z in 1:nrow(starea_grid)) {
  counts <- occ_count(facet="speciesKey", facetLimit=200000,
                      geometry = st_as_text(st_geometry(starea_grid[z,])))
  if (z == 1) {
    counts_list <- counts
  } else {
    counts_list <- bind_rows(counts_list, counts)
    counts_list <- counts_list %>% distinct(speciesKey, .keep_all = T)
  }
  cli_progress_update()
}
cli_progress_done()


# Match both lists
splist_bbox$speciesKey <- as.character(splist_bbox$speciesKey)
match_list <- left_join(counts_list, splist_bbox, by = "speciesKey")

# Filter to remove other ranks and synonyms
# We also remove taxa from kingdoms like Bacteria and Fungi
match_list <- match_list %>%
  filter(taxonRank == "SPECIES") %>%
  filter(taxonomicStatus == "ACCEPTED") %>%
  filter(!kingdom %in% c("Archaea", "Bacteria", "Fungi", "Protozoa", "Viruses")) %>%
  distinct(speciesKey, .keep_all = T)

# See if species is marine
# We need to run in batches due to the way the worms API works...

blocks <- seq(100, nrow(match_list), by = 100)
blocks <- c(blocks, nrow(match_list))

st <- 1
cli_progress_bar(total = length(blocks))
for (z in 1:length(blocks)) {
 
  block_list <- match_list[st:blocks[z],]
  
  worms_res <- try(wm_records_names(block_list$species,
                                    fuzzy = F, marine_only = T))
  
  if (class(worms_res)[1] == "try-error") {
    marine <- FALSE
  } else {
    marine <- TRUE
  }
  
  if (marine) {
    
    worms_res <- lapply(1:length(worms_res), function(x){
      
      wmr <- worms_res[[x]]
      
      if (nrow(wmr) > 0) {
        if ("accepted" %in% wmr$status) {
          wmr <- wmr[wmr$status == "accepted",]
        } else {
          if (!is.na(wmr$valid_AphiaID[1])) {
            wmr <- wm_record(wmr$valid_AphiaID[1])[1,]
          } else {
            wmr <- wmr[1,]
          }
        }
        wmr$gbif_speciesKey <- block_list$speciesKey[x]
        wmr$gbif_taxonKey <- block_list$taxonKey[x]
        wmr$gbif_scientificName <- block_list$scientificName[x]
      }
      return(wmr)
    })
    
    worms_res <- bind_rows(worms_res)
    
    if (!exists("marine_list")) {
      marine_list <- worms_res
    } else {
      marine_list <- bind_rows(marine_list, worms_res)
    }
  }
  
  st <- blocks[z] + 1
  cli_progress_update()
}
cli_progress_done()

# table(marine_list$status)

marine_list$isExtinct[is.na(marine_list$isExtinct)] <- -1

# Filter again to remove those not accepted
marine_list_final <- marine_list %>%
  filter(isExtinct != 1) %>%
  filter(status == "accepted") %>%
  distinct(AphiaID, .keep_all = T) %>%
  filter(rank == "Species")

# Save one with the flags
write.csv(marine_list,
          paste0("data/gbif_splist_flagged_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# Save one without the extincts and with only the unique species
write.csv(marine_list_final,
          paste0("data/gbif_splist_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# END