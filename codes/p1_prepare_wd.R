############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### P1. Prepare working directory ##########################

#> Step 1: Clone the repository to your local machine
#> This will download the scripts and other configuration files

#> Step 2: Download base files from Amazon cloud
#> This will download shapefiles and other accessory files necessary
#download.file("s3://obis.org/mpaeu_sdm", destfile = "amzn_files.zip")
#unzip("amzn_files.zip")
#file.remove("amzn_files.zip")

#> Step 3: Install dependencies
source("requirements.R")

#> Step 4: Check structure
source("check.R")

# END
