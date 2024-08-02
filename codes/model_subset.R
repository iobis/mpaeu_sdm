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
sel_species <- c("Phocoena phocoena", "Arbacia lixula")

# Fit
source("codes/model_fit.R")