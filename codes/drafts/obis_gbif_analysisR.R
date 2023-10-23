#### EMODnet conference - data analysis

# Load packages
library(rgbif)
library(robis)
library(tidyverse)

# Load species list
gbif_list <- read.csv("data/gbif_splist_20230720.csv")
obis_list <- read.csv("data/obis_splist_20230720.csv")

# Ensure that just the same kingdoms are being considered
gbif_list <- gbif_list[gbif_list$kingdom %in% obis_list$kingdom, ]
gbif_list <- gbif_list[!duplicated(gbif_list$gbif_speciesKey),]

# See what is not on OBIS
not_obis <- gbif_list[!gbif_list$valid_AphiaID %in% obis_list$taxonID,]


# Get the counts of points in GBIF to assess rarity
breaks <- seq(1, nrow(not_obis), by = 5)[-1]
breaks[length(breaks)] <- nrow(not_obis)

st <- 1
for (i in breaks) {
  print(st)
  sp_counts <- occ_count(speciesKey = paste(not_obis$gbif_speciesKey[st:i], collapse = ";"), facet=c("speciesKey"), facetLimit=200000)
  #sp_counts$speciesKey <- not_obis$gbif_speciesKey[st:i]
  
  if (st == 1) {
    counts_full <- sp_counts
  } else {
    counts_full <- rbind(counts_full, sp_counts)
  }
  
  st <- i+1
 
}

colnames(counts_full) <- c("gbif_speciesKey", "records_glob_gbif")
counts_full$gbif_speciesKey <- as.integer(counts_full$gbif_speciesKey)

not_obis <- left_join(not_obis, counts_full)

# For now ignore those with NA - synonims, but to check what happens
not_obis <- not_obis[!is.na(not_obis$records_glob_gbif),]

# Get number of records on OBIS
counts_obis <- checklist(taxonid = not_obis$valid_AphiaID)
# Get the counts of points in OBIS
breaks <- seq(1, nrow(not_obis), by = 100)[-1]
breaks[length(breaks)] <- nrow(not_obis)

st <- 1
for (i in breaks) {
  print(st)
  sp_counts <- checklist(taxonid = not_obis$valid_AphiaID[st:i])
  
  if (st == 1) {
    counts_obis <- sp_counts[,c("taxonID", "records")]
  } else {
    counts_obis <- rbind(counts_obis, sp_counts[,c("taxonID", "records")])
  }
  
  st <- i+1
}

colnames(counts_obis) <- c("valid_AphiaID", "records_glob_obis")

not_obis <- left_join(not_obis, counts_obis)

not_obis %>%
  mutate(records_glob_obis = ifelse(is.na(records_glob_obis), 0, records_glob_obis)) %>%
  summarise(rare_gbif = sum(records_glob_gbif <= 10),
            on_obis = sum(records_glob_obis > 0),
            not_on_obis = sum(records_glob_obis == 0),
            on_obis_but_rare = sum(records_glob_obis > 0 & records_glob_obis <= 10))


not_obis %>%
  mutate(records_glob_obis = ifelse(is.na(records_glob_obis), 0, records_glob_obis)) %>%
  filter(records_glob_obis == 0) %>%
  group_by(phylum) %>%
  count()



#### Europe count
sf::sf_use_s2(FALSE)
europe <- sf::st_read("data/shapefiles/mpa_europe_starea_v2.shp")
europe <- sf::st_simplify(sf::st_buffer(europe, dist = 4), dTolerance = 2)

# Get the counts of points in GBIF to assess rarity
breaks <- seq(1, nrow(not_obis), by = 5)[-1]
breaks[length(breaks)] <- nrow(not_obis)

st <- 1
for (i in breaks) {
  print(st)
  sp_counts <- occ_count(speciesKey = paste(not_obis$gbif_speciesKey[st:i], collapse = ";"), facet=c("speciesKey"), facetLimit=200000,
                         geometry = sf::st_as_text(sf::st_as_sfc(europe)))
  #sp_counts$speciesKey <- not_obis$gbif_speciesKey[st:i]
  
  if (st == 1) {
    counts_full <- sp_counts
  } else {
    counts_full <- rbind(counts_full, sp_counts)
  }
  
  st <- i+1
  
}

colnames(counts_full) <- c("gbif_speciesKey", "records_glob_gbif")
counts_full$gbif_speciesKey <- as.integer(counts_full$gbif_speciesKey)

not_obis <- left_join(not_obis, counts_full)



write.csv(not_obis, "obis_gbif_danalisys.csv", row.names = F)
