############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# July of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################# Get habitat info #################################

# Load packages ----
library(obissdm)
library(dplyr)

# Get habitat info ----
species_list <- read.csv(recent_file("data", "all_splist"))

mp_get_ecoinfo(
  species_list$AphiaID,
  outfile = "data/species_ecoinfo.csv",
  overwrite = FALSE,
  return_table = FALSE,
  try_higher = TRUE,
  try_remarks = FALSE,
  show_progress = TRUE
)

# Input with genus info
hab_info <- read.csv("data/species_ecoinfo.csv")
species_list <- species_list %>% 
    select(taxonID, phylum, order, class, family, genus, scientificName)
hab_info <- left_join(hab_info, species_list)

genus <- unique(hab_info$genus)

genus_info <- hab_info %>%
    filter(mode_life != "NOT_FOUND") %>%
    group_by(genus, mode_life) %>%
    count()

hab_info$flag <- 0

# Flags: 
# 0 = found in databases
# 1 = input by genus resolved
# 2 = multiple modes of life in the genus, not possible to resolve
# 3 = genus info not available

for (i in cli::cli_progress_along(genus)) {
    g <- genus[i]
    gen_i <- genus_info %>%
        filter(genus == g) %>%
        mutate(mode_life = case_when(
            grepl("benthic|demersal|bottom", mode_life) ~ "benthic",
            grepl("pelagic$|pelagic_surface|pelagic_mean", mode_life) ~ "pelagic"
        ))
    
    if (nrow(gen_i) > 0) {
        if (length(unique(gen_i$mode_life)) > 1) {
            hab_info$flag[hab_info$genus == g & hab_info$mode_life == "NOT_FOUND"] <- 2
        } else {
            major <- gen_i$mode_life[1]
            hab_info$flag[hab_info$genus == g & hab_info$mode_life == "NOT_FOUND"] <- 1
            hab_info$mode_life[hab_info$genus == g & hab_info$mode_life == "NOT_FOUND"] <- as.character(major)
        }
    } else {
        hab_info$flag[hab_info$genus == g & hab_info$mode_life == "NOT_FOUND"] <- 3
    }
}

table(hab_info$flag)
View(hab_info[hab_info$flag == 2,])

# Try to input by family
# Flags: 
# 0 = found in databases
# 10 = input by family resolved
# 20 = multiple modes of life in the family, not possible to resolve
# 30 = family info not available

hab_info_no_family <- hab_info[is.na(hab_info$family),]
hab_info <- hab_info[!is.na(hab_info$family),]

families <- unique(hab_info$family)

family_info <- hab_info %>%
    filter(mode_life != "NOT_FOUND") %>%
    group_by(family, mode_life) %>%
    count()

for (i in cli::cli_progress_along(families)) {
    f <- families[i]
    fam_i <- family_info %>%
        filter(family == f) %>%
        mutate(mode_life = case_when(
            grepl("benthic|demersal|bottom", mode_life) ~ "benthic",
            grepl("pelagic$|pelagic_surface|pelagic_mean", mode_life) ~ "pelagic"
        ))
    
    if (nrow(fam_i) > 0) {
        if (length(unique(fam_i$mode_life)) > 1) {
            hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] <- 
                hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] + 20
        } else {
            major <- fam_i$mode_life[1]
            hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] <- 
                hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] + 10
            hab_info$mode_life[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] <- as.character(major)
        }
    } else {
        hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] <- 
            hab_info$flag[hab_info$family == f & hab_info$mode_life == "NOT_FOUND"] + 30
    }
}

hab_info <- bind_rows(hab_info, hab_info_no_family)

table(hab_info$flag)

done <- sum(hab_info$flag == 0 | hab_info$flag == 1 | hab_info$flag == 13)
100 - (done * 100)/nrow(hab_info)

# Save new file
write.csv(hab_info[,c("taxonID", "mode_life", "flag")],
    "data/species_ecoinfo.csv", row.names = F)
