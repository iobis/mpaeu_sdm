############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# August of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################## Model a subset of species ###########################

# Settings (remove to use default)
outfolder <- "results"
outacro <- "mpaeu"

# Select a subset to be done
wwf_1 <- read.csv("analysis/wwf_2024/new/wwf_list_082024_matched.csv")
wwf_2 <- read.csv("analysis/wwf_2024/new/wwf_list_082024_birds.csv")

spf <- obissdm::recent_file("data/", "all_sp")
spl <- read.csv(spf)

sel_species <- spl$taxonID[spl$taxonID %in% c(wwf_1$AphiaID, wwf_2$AphiaID)]
sel_species <- sel_species[1:4]
sel_species = 127326

# Fit
source("codes/model_fit.R")