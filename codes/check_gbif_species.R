#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe porject
# s.principe@unesco.org

#### Get list of marine species occurring on the study area on GBIF

# Load packages
library(arrow)
library(sf)
library(terra)
library(tidyverse)
library(worrms)
sf_use_s2(FALSE)

# Load study area
starea <- st_read("data/shapefiles/mpa_europe_starea_v2.shp")
starea <- st_buffer(starea, 0.5)

starea_bbox <- st_bbox(starea)

# Because the study area is quite large, we divide it in blocks to run
# each one separately
grid <- rast(ext(starea_bbox), nrows = 15, ncols = 15)
grid[] <- 1:225
grid <- as.polygons(grid)
grid <- st_as_sf(grid)

starea_grid <- st_intersection(starea, grid)

# Get information from GBIF
# You can both use the AWS full export, or download first a local parquet with all the
# records in the study area. (see get_data_gbif.R)
gbif_snapshot <- "s3://gbif-open-data-eu-central-1/occurrence/2023-06-01/occurrence.parquet"
#gbif_snapshot <- "~/Research/mpa_europe/mpaeu_msdm/data/gbif.parquet"
df <- open_dataset(gbif_snapshot)

for (i in 1:nrow(starea_grid)) {
  cat(i, "\n")
  # Define area
  starea_gr_bbox <- st_bbox(starea_grid[i,])
  
  # Obtain data
  dat <- df %>% 
    select(genus, species, scientificname, decimallongitude, decimallatitude, taxonrank, kingdom) %>%
    #slice_head(n=1) %>%
    filter(
      decimallongitude >= starea_gr_bbox["xmin"] & decimallongitude <= starea_gr_bbox["xmax"],
      decimallatitude >= starea_gr_bbox["ymin"] & decimallatitude <= starea_gr_bbox["ymax"]
    ) %>%
    filter(taxonrank == "SPECIES") %>%
    filter(!kingdom %in% c("Archaea", "Bacteria", "Fungi", "Protozoa")) %>%
    #select(genus, species, scientificname, decimallongitude, decimallatitude) %>%
    collect()
  
  # Convert to SF
  dat_sf <- st_as_sf(dat, coords = c("decimallongitude", "decimallatitude")) 
  st_crs(dat_sf) <- st_crs(starea_grid)
  
  # See if is inside the study area
  inter <- st_contains(starea_grid[i,1], dat_sf, sparse = F)
  
  dat <- dat[inter[1,],]
  
  # Retrieve only unique species
  dat <- dat %>%
    distinct(scientificname, .keep_all = T) %>%
    select(-decimallongitude, -decimallatitude, -taxonrank)
  
  if (i == 1) {
    sp_list <- dat
  } else {
    sp_list <- bind_rows(sp_list, dat)
  }
  rm(dat)
}

# Get only the unique species
sp_list <- sp_list %>%
  distinct(scientificname, .keep_all = T)

# Check which are marine species
breaks <- seq(1, nrow(sp_list), by = 50)
breaks <- c(breaks, nrow(sp_list))
breaks <- breaks[-1]
st <- 1

# Get one example of result to use when nothing is returned
ex <- wm_records_names("Acanthurus chirurgus")
ex <- ex[[1]]
ex[,] <- NA

# Run for all
cli::cli_progress_bar("Checking species", total = length(breaks))
for (i in 1:length(breaks)) {
  recs <- try(wm_records_names(sp_list$species[st:breaks[i]]))
  if (class(recs) != "try-error") {
    recs <- lapply(1:length(recs), function(z){
      x <- recs[[z]]
      if (nrow(x) < 1) {
        x <- ex
      } else {
        if (nrow(x) > 1) {
          x <- x[x$status != "unaccepted",]
          if (nrow(x) > 1) {
            if ("accepted" %in% x$status) {
              x <- x[x$status == "accepted",]
            }
          }
          if (nrow(x) > 1) {
            auth <- sp_list$scientificname[st:breaks[i]][z]
            auth <- str_trim(gsub("\\(|\\)", "", auth))
            full_names <- mutate(x, f_nam = paste(scientificname, authority))
            ldist <- adist(auth, full_names$f_nam)[1,]
            x <- x[which.min(ldist),]
            warning("Multiple matches for ", auth, ".\nMatched with: ", x$scientificname, " ", x$authority)
          }
          if (nrow(x) < 1) {
            x <- ex
          }
        } else {
          x
        }
      }
      return(x)
    })
    recs <- bind_rows(recs)
    recs$gbif_species <- sp_list$species[st:breaks[i]]
    recs$gbif_names <- sp_list$scientificname[st:breaks[i]]
    if (i == 1) {
      full_recs <- recs
    } else {
      full_recs <- bind_rows(full_recs, recs)
    }
  } else {
    recs <- ex
    recs <- tibble::add_row(recs, AphiaID = rep(NA, length(st:breaks[i])-1))
    recs$gbif_species <- sp_list$species[st:breaks[i]]
    recs$gbif_names <- sp_list$scientificname[st:breaks[i]]
    if (i == 1) {
      full_recs <- recs
    } else {
      full_recs <- bind_rows(full_recs, recs)
    }
  }
  st <- breaks[i]+1
  cli::cli_progress_update()
}
cli::cli_progress_done()

# Remove those that are not marine
final_list <- full_recs[!is.na(full_recs$AphiaID),]
final_list$isExtinct <- ifelse(is.na(final_list$isExtinct),
                               -1, final_list$isExtinct)

# Save one with the flags
write.csv(final_list,
          paste0("data/gbif_splist_flagged_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# Save one without the extincts and with only the unique species
final_list <- final_list[final_list$isExtinct != 1, ]
final_list <- final_list %>%
  distinct(valid_name, .keep_all = T)

write.csv(final_list,
          paste0("data/gbif_splist_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

# END