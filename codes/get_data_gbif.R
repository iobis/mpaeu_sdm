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
  splist$scientificName, splist$scientificNameAuthorship
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
check_nm <- name_backbone_checklist(non_match$scientificName,
                                    phylum = non_match$phylum,
                                    class = non_match$class,
                                    family = non_match$family,
                                    genus = non_match$genus,
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



# Assign checking using synonim
check_syn_nm <- lapply(non_match$taxonID, function(x){
  res <- try(wm_synonyms(x))
  nv <- NA
  if (class(res)[1] != "try-error") {
    if (any(res$status %in% c("alternate representation", "superseded combination"))) {
      syn <- res$scientificname[res$status %in%  c("alternate representation", "superseded combination")][1]
      cn <- name_backbone(syn, strict = T)
      if (cn$matchType != "NONE") {
        nv <- cn$usageKey[1]
      }
    }
  }
  Sys.sleep(2)
  return(nv)
})

check_syn_nm <- unlist(check_syn_nm)

non_match$usageKey <- check_syn_nm

# Match synonims
splist_match <- left_join(splist_match, non_match[,c("taxonID", "usageKey")],
                          by = "taxonID")
splist_match$usageKey <- splist_match$usageKey.x
splist_match$usageKey[!is.na(splist_match$usageKey.y)] <- splist_match$usageKey.y[!is.na(splist_match$usageKey.y)]

splist_match <- select(splist_match, -usageKey.x, -usageKey.y, -full_name)

splist_match <- splist_match[!is.na(splist_match$usageKey),]

# Get data from GBIF and save
mp_get_gbif(sci_names = splist_match$usageKey)

# END