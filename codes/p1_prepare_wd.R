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
# The easiest way is downloading the AWS CLI program: https://aws.amazon.com/cli/
# Then use this on the command line:
# aws s3 cp --recursive s3://obis-maps/sdm/source/model=mpaeu/ . --no-sign-request

#> Step 3: Install dependencies
source("requirements.R")

#> Step 4: Check structure
source("check.R")

# END
