############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### P2. Download necessary data ############################

#> Step 1: Obtain species list
#> Two ways of downloading species lists are available. Since June 2024 we are
#> applying the version 2, that uses the OBIS gridded product
# For the previous version, see `get_species_lists.R`
source("codes/get_species_lists_grid.R")

#> Step 2: Download species data
#> This will download species data from OBIS and GBIF
# By default, for OBIS it will download the full export. 
# An alternative option is to download just the subset. See `get_species_data.R`
source("codes/get_species_data_full.R")

#> Step 3: Download environmental data
#> This will download environmental data from Bio-ORACLE and prepare other files
source("codes/get_env_data.R")
source("codes/get_add_env_data.R")

#> Step 5: Get habitat information
#> Obtain habitat information from WoRMS, SeaLifeBase and FishBase
source("codes/get_habitat_information.R")

# END