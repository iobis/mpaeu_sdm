############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
###################### EMODNet Open Conference 2023 ###########################
# Generate figure for the poster

library(treemap)
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

# Put together in a data frame
dat <- data.frame(
  where = as.factor(c("GBIF", "OBIS", "Shared")),
  n = c(nrow(not_obis), nrow(not_gbif), nrow(both_lists))
)

# Plot
treemap(dat,
        # data
        index="where",
        vSize="n",
        type="index",
        # Main
        title="",
        palette=c("#94d2bd", "#2a9d8f", "#005f73"),
        # Borders:
        border.col=c("white"),             
        border.lwds=1,                       
        # Labels
        # fontsize.labels=1,
        fontcolor.labels="white",
        fontface.labels=1,            
        # bg.labels=c("transparent"),              
        align.labels=c("left", "top")                                  
        # overlap.labels=0.5,
        # inflate.labels=T
        )

# Save as PNG
png("treemap_emodnetposter.png", width = 600, height = 350)
treemap(dat,
        # data
        index="where",
        vSize="n",
        type="index",
        # Main
        title="",
        palette=c("#94d2bd", "#2a9d8f", "#005f73"),
        # Borders:
        border.col=c("white"),             
        border.lwds=1,                       
        # Labels
        # fontsize.labels=1,
        fontcolor.labels="white",
        fontface.labels=1,            
        # bg.labels=c("transparent"),              
        align.labels=c("left", "top")                                  
        # overlap.labels=0.5,
        # inflate.labels=T
)
dev.off()
