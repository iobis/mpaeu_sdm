# Get information from Worms about traits
library(tidyverse)
library(worrms)

splist_gbif <- read.csv("data/gbif_splist_20230720.csv")
splist_obis <- read.csv("data/obis_splist_20230720.csv")

splist_gbif <- splist_gbif[,c("AphiaID", "scientificname")]
splist_obis <- splist_obis[,c("taxonID", "scientificName")]

colnames(splist_gbif) <- colnames(splist_obis) <- c("AphiaID", "scientificName")

splist <- rbind(splist_gbif, splist_obis)

splist <- distinct(splist, AphiaID, .keep_all = T)

###
teste <- wm_attr_data(splist$AphiaID[1], include_inherited = T)
