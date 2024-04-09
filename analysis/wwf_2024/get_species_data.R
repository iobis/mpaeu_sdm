# Prepare WWF list of species

# Load WWF list
wwf_list <- readxl::read_excel("~/Downloads/Species List.xlsx")
colnames(wwf_list) <- c(
  "vernacularName",
  "scientificNameWWF",
  "rank",
  "listedIn",
  "annex",
  "note"
)

head(wwf_list)

# Match with WoRMS
new_names <- obistools::match_taxa(wwf_list$scientificNameWWF)
# For Patella ferruginea we go with the accepted name 140679

# Bind new names
new_wwf_list <- bind_cols(new_names, wwf_list)

new_wwf_list <- new_wwf_list %>%
  mutate(AphiaID = acceptedNameUsageID) %>%
  relocate(AphiaID, scientificName, scientificNameWWF, vernacularName)

# Correct the Pusa record with new name
pusa <- worrms::wm_record(159021)
colnames(pusa)[3] <- "scientificName"
pusa$scientificNameWWF <- "Phoca hispida"

new_wwf_list_b <- bind_rows(new_wwf_list %>% 
                              filter(!stringr::str_detect(scientificNameWWF, "Phoca hispida")) %>%
                              mutate(AphiaID = as.numeric(AphiaID)),
                            pusa)
new_wwf_list_b <- new_wwf_list_b[,colnames(new_wwf_list)]

# Get Alosa records
alosa <- worrms::wm_children(125715)
alosa <- alosa %>%
  filter(status == "accepted") %>% 
  filter(rank == "Species")

# Get Cetaceans records
cetaceans <- worrms::wm_children(2688)
cetaceans <- cetaceans %>% filter(status == "accepted")
# Now Families
cetaceans <- worrms::wm_children_(cetaceans$AphiaID)
cetaceans <- cetaceans %>% filter(status == "accepted")
# Now Genus
cetaceans <- worrms::wm_children_(cetaceans$AphiaID)
cetaceans <- cetaceans %>% filter(status == "accepted")
# Now species
cetaceans <- worrms::wm_children_(cetaceans$AphiaID)
cetaceans_acc <- cetaceans %>%
  filter(status == "accepted") %>%
  filter(rank == "Species")

# Get Acipenseridae
acipenseridae <- worrms::wm_children(worrms::wm_name2id("Acipenseridae"))
acipenseridae <- acipenseridae %>%
  filter(status == "accepted")
acipenseridae <- worrms::wm_children_(acipenseridae$AphiaID)
acipenseridae <- acipenseridae %>%
  filter(status == "accepted") %>%
  filter(rank == "Species")

# Get Rhinobatidae
rhinobatidae <- worrms::wm_children(worrms::wm_name2id("Rhinobatidae"))
rhinobatidae <- rhinobatidae %>%
  filter(status == "accepted")
rhinobatidae <- worrms::wm_children_(rhinobatidae$AphiaID)
rhinobatidae <- rhinobatidae %>%
  filter(status == "accepted") %>%
  filter(rank == "Species")

# Correct column names
cetaceans <- cetaceans %>%
  rename(scientificName = scientificname)
rhinobatidae <- rhinobatidae %>%
  rename(scientificName = scientificname)
acipenseridae <- acipenseridae %>%
  rename(scientificName = scientificname)
alosa <- alosa %>%
  rename(scientificName = scientificname)

# Bind new names
new_wwf_list_c <- new_wwf_list_b %>%
  mutate(status = "accepted") %>%
  filter(scientificNameWWF != "Cetaceans") %>%
  filter(scientificName != "Alosa") %>%
  filter(scientificName != "Rhinobatidae") %>%
  filter(scientificNameWWF != "Acipenseridae") %>%
  bind_rows(cetaceans, alosa, acipenseridae, rhinobatidae) %>%
  filter(status == "accepted") %>%
  filter(rank == "Species")

# Add GBIF information
gbif_names <- lapply(new_wwf_list_c$scientificName, rgbif::name_backbone)
gbif_names <- bind_rows(gbif_names)
gbif_names <- gbif_names %>%
  select(gbif_speciesKey = usageKey, gbif_scientificName = scientificName, gbif_matchType = matchType)

new_wwf_list_d <- new_wwf_list_c %>% bind_cols(gbif_names)

# Remove duplicates
new_wwf_list_final <- new_wwf_list_d %>%
  distinct(AphiaID, .keep_all = T)

new_wwf_list_final$taxonID <- new_wwf_list_final$AphiaID

# Fill details of those missing
missing_det <- worrms::wm_record_(id = new_wwf_list_final$AphiaID[is.na(new_wwf_list_final$phylum)])
missing_det <- bind_rows(missing_det)
missing_det <- missing_det[,c("kingdom", "phylum", "class", "order", "family")]

new_wwf_list_final[is.na(new_wwf_list_final$phylum), c("kingdom", "phylum", "class", "order", "family")] <- missing_det

write.csv(new_wwf_list_final, "analysis/wwf_2024/wwf_list.csv", row.names = F)

# Download data not available on our previous download
all_sp <- read.csv("data/all_splist_20240319.csv")

not_available <- new_wwf_list_final[!new_wwf_list_final$AphiaID %in% all_sp$taxonID,]

setwd("../mpaeu_shared/")
obissdm::mp_get_gbif(sci_names = not_available$gbif_speciesKey)
setwd("../mpaeu_sdm/")

# And then do the standardization
for (i in 1:length(unique(not_available$AphiaID))) {
  if (!file.exists(paste0("data/species/", "key=", unique(not_available$AphiaID)[i], ".parquet"))) {
    cat(unique(not_available$AphiaID)[i], "\n")
    mp_standardize(species = unique(not_available$AphiaID)[i],
                   species_list = "analysis/wwf_2024/wwf_list.csv",
                   sdm_base = rast("data/env/current/thetao_baseline_depthsurf_mean.tif"),
                   species_folder = "../mpaeu_shared/")
  } else {
    cat(unique(not_available$AphiaID)[i], "done \n")
  }
}

# END