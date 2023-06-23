#### MPA EUROPE PROJECT ####
# June 2023
# OBIS contribution to the MPA Europe porject
# s.principe@unesco.org

### Get the data from GBIF
# Load packages
library(rgbif)
library(worrms)
library(obissdm)

# Load species list
splist <- read.csv("data/gbif_splist_20230623.csv")

splist$full_name <- paste(
  splist$scientificname, splist$authority
)
  
# Check GBIF taxonomic backbone
check_names <- name_backbone_checklist(splist$full_name,
                                       strict = T)
check_names$taxonID <- splist$taxonID

# See if any name was not matched
sum(is.na(check_names$usageKey))

non_match <- splist[is.na(check_names$usageKey),]

# Manually correct those
# Check using only the species name
check_nm <- name_backbone_checklist(non_match$scientificname,
                                    strict = T)
check_nm$taxonID <- non_match$taxonID

# Bind those that are already ok
splist_match <- left_join(splist, check_names[,c("taxonID", "usageKey")],
                          by = "taxonID")
splist_match <- left_join(splist_match, check_nm[,c("taxonID", "usageKey")],
                          by = "taxonID")
splist_match$usageKey <- splist_match$usageKey.x
splist_match$usageKey[!is.na(splist_match$usageKey.y)] <- splist_match$usageKey.y[!is.na(splist_match$usageKey.y)]

splist_match <- select(splist_match, -usageKey.x, -usageKey.y)

non_match <- check_nm[is.na(check_nm$usageKey),]

# Get data from GBIF and save
mp_get_gbif(sci_names = splist_match$usageKey)

# END