library(obissdm)
library(duckdb)
library(sf)
library(glue)
library(arrow)
library(dplyr)
# For the optional part:
library(terra)


# Prepare settings and objects -----
con <- dbConnect(duckdb())
dbSendQuery(con, "install spatial; load spatial;")

study_area <- st_read("data/shapefiles/mpa_europe_starea_v3.gpkg")

study_wkt <- st_as_text(st_as_sfc(study_area))

species_list <- read.csv(recent_file("data", "all_splist"))


# Check raw data -----
obis_src <- "data/raw/obis_20240625.parquet"
gbif_src <- "data/raw/gbif_20240726/*/*.parquet"

obis_results <- dbGetQuery(con, glue(
"
select 
    AphiaID,
    scientificName,
    sum(case when ST_Intersects(ST_Point(decimalLongitude, decimalLatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_intersects,
    sum(case when not ST_Intersects(ST_Point(decimalLongitude, decimalLatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_non_intersects
from read_parquet('{obis_src}')
where basisOfRecord not in ('FossilSpecimen', 'FOSSIL_SPECIMEN', 'LivingSpecimen', 'LIVING_SPECIMEN') and date_year >= 1950
group by AphiaID, scientificName
"
))

gbif_results <- dbGetQuery(con, glue(
"
select 
    specieskey as gbif_speciesKey,
    species as gbif_scientificName,
    sum(case when ST_Intersects(ST_Point(decimallongitude, decimallatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_intersects,
    sum(case when not ST_Intersects(ST_Point(decimallongitude, decimallatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_non_intersects
from read_parquet('{gbif_src}')
where basisofrecord not in ('FossilSpecimen', 'FOSSIL_SPECIMEN', 'LivingSpecimen', 'LIVING_SPECIMEN') and year >= 1950
group by specieskey, species
"
))

# Check processed data ------
proc_src <- "data/species/*.parquet"

proc_results <- dbGetQuery(con, glue(
    "
select 
    taxonID,
    data_type,
    sum(case when ST_Intersects(ST_Point(decimallongitude, decimallatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_intersects,
    sum(case when not ST_Intersects(ST_Point(decimallongitude, decimallatitude), ST_GeomFromText('{study_wkt}')) then 1 else 0 end) as n_non_intersects
from read_parquet('{proc_src}')
group by taxonID, data_type
    "
))

proc_results <- proc_results |>
    mutate(n_total = n_intersects + n_non_intersects) |>
    rename(on_area = n_intersects, out_area = n_non_intersects, total = n_total) |>
    tidyr::pivot_wider(names_from = data_type, values_from = c(on_area, out_area, total))

dbDisconnect(con)


# Join results -----
obis_results <- obis_results |>
    select(taxonID = AphiaID,
           obis_on_area = n_intersects,
           obis_out_area = n_non_intersects)

species_list_ed <- species_list |>
    left_join(obis_results, by = "taxonID")

gbif_results <- gbif_results |>
    select(gbif_speciesKey,
           gbif_on_area = n_intersects,
           gbif_out_area = n_non_intersects)

species_list_ed <- species_list_ed |>
    left_join(gbif_results, by = "gbif_speciesKey")

species_list_ed <- species_list_ed |>
    left_join(proc_results, by = "taxonID")


# Optional: check area that overlay with study area on processed rasters -----
# This is generated through codes/post_prep_layers_zonation.R
proc_rasters <- list.files("results/processed_layers", full.names = T)

proc_rasters <- proc_rasters[grepl("scen=current_th=p10_type=const", proc_rasters)]

pb <- progress::progress_bar$new(total = length(proc_rasters))
covered_areas <- data.frame(
    taxonID = as.integer(gsub("taxonid=", "", gsub("_.*", "", basename(proc_rasters)))),
    present_area_perc = NA
)
for (i in seq_along(proc_rasters)) {
    pb$tick()
    layer <- rast(proc_rasters[i])
    layer <- classify(layer, matrix(c(0, Inf, 1), nrow = 1))
    layer_a <- expanse(layer, unit = "km", byValue = T)
    if (length(layer_a$value) > 1) {
        pres <- layer_a$area[layer_a$value == 1]
        pres <- (pres * 100) / sum(layer_a$area)
    } else if (layer_a$value == 0) {
        pres <- 0
    } else if (layer_a$value == 1) {
        pres <- 100
    } else {
        stop("Check layer.")
    }
    covered_areas$present_area_perc[i] <- round(pres, 1)
}

species_list_ed <- species_list_ed |>
    left_join(covered_areas, by = "taxonID")


# Optional: add SDM status -------
sdm_status <- read.table("internal/sdm_status_20250925.txt", header = T)
table(sdm_status$sdm_status)
sdm_status$sdm_status[sdm_status$sdm_status == "all_failed"] <- "failed"

species_list_ed <- species_list_ed |>
    left_join(sdm_status, by = "taxonID")
species_list_ed$sdm_status[is.na(species_list_ed$sdm_status)] <- "not_available"


# Save object -----
write.csv(species_list_ed,
          paste0("internal/rarity_check_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = FALSE)


# Example analysis -----
species_table <- read.csv("internal/rarity_check_20250925.csv")

# Number of species that had records in the study area after data cleaning:
species_table |>
    mutate(on_area_fit_points = ifelse(is.na(on_area_fit_points), 0, on_area_fit_points)) |>
    summarise(
        with_data = sum(on_area_fit_points > 0),
        no_data = sum(on_area_fit_points == 0)
    )

# From those with data, the number of species that had at least 15 records 
# (enough for ESM) or at least 30 records (for standard SDM):
species_table |>
    filter(!is.na(on_area_fit_points)) |>
    filter(on_area_fit_points > 0) |>
    summarise(
        low_data = sum(on_area_fit_points < 15),
        data_for_esm = sum(on_area_fit_points >= 15 & on_area_fit_points < 30),
        data_for_standard = sum(on_area_fit_points >= 30)
    )

# Get which species were modelled but have small percentage of suitable area covering
# the study area
species_low_area <- species_table |>
    filter(sdm_status == "succeeded") |>
    filter(present_area_perc < 5) |> # Trying 5% as a threshold, but can be anything....
    select(taxonID, scientificName, percentage_suitable_on_studyarea = present_area_perc)

head(species_low_area)
nrow(species_low_area)
