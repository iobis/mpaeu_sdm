fao_areas_fb <- rfishbase::fb_tbl(
  "faoareas",
  server = c("fishbase"),
  version = "latest",
  collect = TRUE
)

fao_areas_sb <- rfishbase::fb_tbl(
  "faoareas",
  server = c("sealifebase"),
  version = "latest",
  collect = TRUE
)

fao_areas_fb$DateEntered <- as.character(fao_areas_fb$DateEntered)
fao_areas_fb$DateModified <- as.character(fao_areas_fb$DateModified)
fao_areas_fb$DateChecked <- as.character(fao_areas_fb$DateChecked)
fao_areas_fb$TS <- as.character(fao_areas_fb$TS)

fao_area <- bind_rows(fao_areas_fb, fao_areas_sb)

arrow::write_parquet(fao_area, "data/fao_areas.parquet")
###


obis_list <- read.csv("data/obis_splist_20230720.csv")
gbif_list <- read.csv("data/gbif_splist_20231023.csv")

obis_list <- obis_list %>%
  mutate(origin = "obis")

gbif_list <- gbif_list %>%
  mutate(species = scientificname) %>%
  rename(taxonID = AphiaID) %>%
  mutate(origin = "gbif")

all_list <- bind_rows(obis_list, gbif_list)

all_list <- all_list %>%
  distinct(taxonID, .keep_all = T) %>% 
  filter(kingdom != "Viruses" & kingdom != "Protozoa")

write.csv(all_list, paste0("data/all_splist_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)


###

teste_sp <- 124316

mp_standardize(species = 124316,
               sdm_base = rast("data/env/current/thetao_baseline_depthsurf_mean.tif"),
               species_list = "data/all_splist_20240319.csv",
               species_folder = "../mpaeu_shared/")
