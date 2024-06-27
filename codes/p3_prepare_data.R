############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################### P3. Prepare data for use on models #########################

#> Step 1: Standardize datasets
#> Perform automated QC steps on the data and convert to format
source("codes/prepare_dataqc.R")

#> Step 2: Test environmental data for colinearity
# Note that this should be done manually. Running the source code here
# will just apply the same steps used in the original project
source("codes/prepare_envdata.R")

# END