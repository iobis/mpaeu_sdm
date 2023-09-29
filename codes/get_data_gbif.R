#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe porject
# s.principe@unesco.org

### Get the data from GBIF
# Load packages
library(obissdm)

# Load species list
splist <- read.csv("data/gbif_splist_20230720.csv")

# Because we already have the GBIF taxonKey we can supply this to the function

# Get data from GBIF and save
mp_get_gbif(sci_names = splist$gbif_taxonKey)

# END