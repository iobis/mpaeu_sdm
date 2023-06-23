#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe porject
# s.principe@unesco.org

#### Get list of marine species occurring on the study area on GBIF

# Load packages
library(sf)
library(terra)
library(tidyverse)
library(robis)
library(worrms)
sf_use_s2(FALSE)

# Load study area
starea <- st_read("data/shapefiles/mpa_europe_starea_v2.shp")
starea <- st_buffer(starea, 0.5)

starea_bbox <- st_bbox(starea)

# Because the study area is quite large, we divide it in blocks to run
# each one separately
grid <- rast(ext(starea_bbox), nrows = 25, ncols = 25)
grid[] <- 1:625
grid <- as.polygons(grid)
grid <- st_as_sf(grid)

starea_grid <- st_intersection(starea, grid)

# Get the checklist of species for each area
for (i in 1:nrow(starea_grid)) {
  cat(i)
  dat <- checklist(geometry = st_as_text(st_geometry(starea_grid[i,])))
  if (i == 1) {
    sp_list <- dat
  } else {
    sp_list <- bind_rows(sp_list, dat)
  }
  rm(dat)
}

# Filter according to project
final_list <- sp_list %>%
  filter(taxonRank == "Species") %>%
  filter(!kingdom %in% c("Archaea", "Bacteria", "Fungi", "Protozoa")) %>%
  filter(taxonomicStatus == "accepted") %>% 
  filter(!is.na(is_marine) | !is.na(is_brackish) | !is.na(is_terrestrial) | !is.na(is_freshwater)) %>%
  rowwise() %>%
  mutate(btm_na = ifelse(sum(is.na(c(is_brackish, is_terrestrial, is_marine))) < 3, 1, NA)) %>%
  mutate(is_freshwater = ifelse(is.na(is_freshwater), FALSE, is_freshwater)) %>%
  filter(!is_freshwater & !is.na(btm_na)) %>%
  filter(!is.na(phylum)) %>%
  distinct(scientificName, .keep_all = T)

# Check which ones are extant species by getting the full information from WoRMS
# We need to do in chunks, otherwise the server returns an error
breaks <- seq(1, nrow(final_list), by = 49)
breaks <- c(breaks, nrow(final_list))
breaks <- breaks[-1]
st <- 1

cli::cli_progress_bar("Checking species", total = length(breaks))
for (i in 1:length(breaks)) {
  info <- wm_record(final_list$taxonID[st:breaks[i]])
  if (nrow(info) < length(final_list$taxonID[st:breaks[i]])) {
    stop("Not equal number retrieved")
  }
  if (i == 1) {
    is_ext <- info[,c("AphiaID", "isExtinct")]
  } else {
    is_ext <- bind_rows(is_ext, info[,c("AphiaID", "isExtinct")])
  }
  st <- breaks[i]+1
  cli::cli_progress_update()
}
cli::cli_progress_done()

# Add information to the main table
colnames(is_ext)[1] <- "taxonID"

# Save a flagged version
final_list <- left_join(final_list, is_ext, by = "taxonID")
final_list$isExtinct[is.na(final_list$isExtinct)] <- -1

write.csv(final_list,
          paste0("data/obis_splist_flagged_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# Save one without the extincts
final_list <- final_list[final_list$isExtinct != 1, ]

write.csv(final_list,
          paste0("data/obis_splist_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# END