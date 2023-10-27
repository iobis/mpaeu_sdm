############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2023
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############### Generate standardized species data snapshot ####################
# Not all data produced by this project cna be easily shared by GitHub due to size
# limitations. Specially, the species information is huge, because we have files
# for OBIS and GBIF data, a file for unduplicated data and finally a file for
# all standardized data.
# 
# To facilitate the use of this data, we provide one snapshot of the most up-to-date
# standardized species data (produced using function obissdm::mp_standardize())
# This snapshot, saved as a parquet file, contains only data for species with more
# than 30 occurrence records

# Load packages ----
library(arrow)
library(dplyr)

# Prepare data ----
# Get first time to read schema
first <- read_parquet("data/species/key=100805/date=20231022/ftype=stdpts/spdata0.parquet")
first_sch <- schema(arrow_table(first))
first_sch$ftype = field("ftype", string())

# Open dataset based on schema
ds <- open_dataset("data/species/", schema = first_sch)

# Open dataset and filter
final <- ds %>%
  filter(ftype == "stdpts") %>%
  collect()

# Get counts to filter those with less than 30
counts <- final %>%
  group_by(species) %>%
  count()

final <- final %>%
  filter(species %in% counts$species[counts$n >= 30])

# Save final file ----
# Ensure folder exists
fs::dir_create("snapshot")
# If any other file exists on the folder, remove it
if (length(list.files("snapshot")) > 0) {
  file.remove(list.files("snapshot", full.names = T))
}
write_parquet(final, paste0("snapshot/std_records_", format(Sys.Date(), "%Y%m%d"), ".parquet"))

# Save also final species list file
sp_list <- final %>%
  distinct(taxonID)

# Load other information about species
full_info <- read.csv("data/all_splist_20231017.csv")
colnames(full_info)[1] <- "taxonID"

sp_list <- left_join(sp_list, full_info)

write_parquet(sp_list, paste0("snapshot/std_splist_", format(Sys.Date(), "%Y%m%d"), ".parquet"))
              