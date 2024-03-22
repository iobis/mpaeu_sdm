
sel_species <- list.files("data/species/")
sel_species <- gsub("key=", "", sel_species)
sel_species <- gsub("\\.parquet", "", sel_species)
sel_species <- as.numeric(sel_species)

source("codes/model_fit.R")