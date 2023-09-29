par(mfrow = c(1,1))
plot(original_suit_102);points(sp_data_102, pch = 20, cex = .5)
par(mfrow = c(1,2))
plot(original_suit_102);plot(pred.max_102, main = "Maxnet")
plot(original_suit_102);plot(pred.brt_102, main = "BRT")
plot(original_suit_102);plot(pred.rf_102, main = "RF")
plot(original_suit_102);plot(pred.lasso_102, main = "LASSO")



par(mfrow = c(1,1))
plot(original_suit_101);points(sp_data_101, pch = 20, cex = .5)
par(mfrow = c(1,2))
plot(original_suit_101);plot(pred.max_101, main = "Maxnet")
plot(original_suit_101);plot(pred.brt_101, main = "BRT")
plot(original_suit_101);plot(pred.rf_101, main = "RF")
plot(original_suit_101);plot(pred.lasso_101, main = "LASSO")

teste.max

plot(resp_curves(teste.rf, env))


library(dplyr)

tabela_gbif <- read.csv("../mpaeu_sdm/data/gbif_splist_20230720.csv")
tabela_gbif <- tabela_gbif %>%
  filter(!valid_AphiaID %in% c(816335, 1486516))

tabela_obis <- read.csv("../mpaeu_sdm/data/obis_splist_20230720.csv")

unique_gbif <- unique(tabela_gbif$valid_AphiaID)
unique_obis <- unique(tabela_obis$taxonID)

shared_sp <- sum(unique_gbif %in% unique_obis)
only_gbif <- sum(!unique_gbif %in% unique_obis)
only_obis <- sum(!unique_obis %in% unique_gbif)

sum(shared_sp, only_gbif, only_obis) == length(unique(c(unique_gbif, unique_obis)))


species_gbif <- tabela_gbif %>%
  select(aphiaid = valid_AphiaID, species = scientificname, phylum) %>%
  mutate(origin = "GBIF")

species_obis <- tabela_obis %>%
  select(aphiaid = taxonID, species = scientificName, phylum) %>%
  mutate(origin = "OBIS")

#all_species <- bind_rows(species_gbif, species_obis)
all_species <- full_join(species_gbif, species_obis, by = "aphiaid")

all_species <- all_species %>%
  rowwise() %>%
  mutate(where = ifelse(
    !is.na(origin.x) & !is.na(origin.y), "BOTH",
    ifelse(is.na(origin.y), "GBIF", "OBIS")
  )) %>%
  mutate(phylumb = ifelse(
    is.na(phylum.x), phylum.y, phylum.x
  ))

group_count <- all_species %>%
  group_by(phylumb, where) %>%
  count()

sum(group_count$n)

write.csv(group_count, "~/Research/records_byphylum_bydatabase.csv", row.names = F)
