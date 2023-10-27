#### EMODnet conference - data analysis

# Load packages
library(rgbif)
library(robis)
library(tidyverse)

# Load species list
gbif_list <- read.csv("data/gbif_splist_20231023.csv") # Most recent, with improvement on code
obis_list <- read.csv("data/obis_splist_20230720.csv")

# Ensure that just the same kingdoms are being considered
gbif_list <- gbif_list[gbif_list$kingdom %in% obis_list$kingdom, ]
gbif_list <- gbif_list[!duplicated(gbif_list$gbif_speciesKey),]

# See what is not on OBIS
not_obis <- gbif_list[!gbif_list$valid_AphiaID %in% obis_list$taxonID,]
nrow(not_obis)
# See what is not on GBIF
not_gbif <- obis_list[!obis_list$taxonID %in% gbif_list$valid_AphiaID,] 
nrow(not_gbif)
# See what is shared
both_lists <- obis_list[obis_list$taxonID %in% gbif_list$valid_AphiaID,]
nrow(both_lists)


# Get the counts of points in GBIF globally to assess rarity
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

# For now ignore those with NA - synonims for GBIF, but not for worms...
not_obis <- not_obis[!is.na(not_obis$records_glob_gbif),]


# Get the counts of points in OBIS globally
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


# We can then see which records are missing on OBIS just in Europe or globally
# as well as seeing which species are also rare on GBIF.
not_obis_summary <- not_obis %>%
  mutate(records_glob_obis = ifelse(is.na(records_glob_obis), 0, records_glob_obis)) %>%
  summarise(rare_gbif = sum(records_glob_gbif <= 10),
            on_obis = sum(records_glob_obis > 0),
            not_on_obis = sum(records_glob_obis == 0),
            on_obis_but_rare = sum(records_glob_obis > 0 & records_glob_obis <= 10))

not_obis_summary

not_obis_phylum <- not_obis %>%
  mutate(records_glob_obis = ifelse(is.na(records_glob_obis), 0, records_glob_obis)) %>%
  filter(records_glob_obis == 0) %>%
  group_by(phylum) %>%
  count()

not_obis_phylum <- bind_rows(not_obis_phylum,
                             data.frame(
                               phylum = "Others",
                               n = sum(not_obis_phylum$n[not_obis_phylum$n < 2])
                             ))
not_obis_phylum <- not_obis_phylum %>% filter(n >= 2)



ggplot(not_obis_phylum) +
  geom_bar(aes(x = reorder(phylum, -n), y = n), stat = "identity", fill = "#006dd7") +
  xlab("Phylum") + ylab("Number of species missing on OBIS") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.x = element_blank())
ggsave("records_by_phylum.png", width = 8, height = 6)


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

colnames(counts_full) <- c("gbif_speciesKey", "records_europe_simp_gbif")
counts_full$gbif_speciesKey <- as.integer(counts_full$gbif_speciesKey)

not_obis <- left_join(not_obis, counts_full)



write.csv(not_obis, "obis_gbif_danalisys.csv", row.names = F)
