#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe porject
# s.principe@unesco.org

### Get the data from OBIS
# Load packages
library(robis)
library(obissdm)

# Load species list
splist <- read.csv("data/obis_splist_20230622.csv")

# Get data from OBIS and save
mp_get_obis(sci_names = splist$taxonID)

# END