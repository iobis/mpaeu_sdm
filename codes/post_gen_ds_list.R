############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# February of 2025
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################# Generate list of datasets used in the work ###################

# Load packages/settings -----
library(arrow)
library(dplyr)
version <- "1.0"
cur_date <- format(Sys.Date(), "%Y-%m-%d")

obis_ds <- "data/raw/obis_20240625.parquet"
gbif_ds <- "data/raw/gbif_20240726/"

sp_list <- read.csv("data/all_splist_20240724.csv")


# List used datasets
obis_ds <- open_dataset(obis_ds)

obis_list <- obis_ds |>
    filter(AphiaID %in% sp_list$taxonID) |>
    group_by(AphiaID, dataset_id) |>
    count() |> collect()

gbif_ds <- open_dataset(gbif_ds)

gbif_keys <- sp_list$gbif_speciesKey[!is.na(sp_list$gbif_speciesKey)]

gbif_list <- gbif_ds |>
    filter(specieskey %in% gbif_keys) |>
    group_by(specieskey, datasetkey) |>
    count() |> collect()

head(obis_list)
head(gbif_list)

obis_list$source <- "obis"
gbif_list$source <- "gbif"

colnames(gbif_list) <- colnames(obis_list) <- c("key", "dataset_id", "records", "source")

pb <- progress::progress_bar$new(total = length(unique(gbif_list$dataset_id)))
gbif_ds_context <- lapply(unique(gbif_list$dataset_id), function(di) {
    pb$tick()
    dc <- rgbif::dataset_get(di)
    dc |> 
        select(dataset_id = "key", any_of(c("doi", "title", "citation", "description")))  |> 
        (\(x) {
            if ("citation" %in% colnames(x)) {
                mutate(x, citation = paste(unlist(citation)[c("text", "identifier")], collapse = "|"))
            } else {
                x
            }
        })()
})

obis_ds_context <- robis::dataset(datasetid = unique(obis_list$dataset_id)[1:100])

obis_ids <- split(unique(obis_list$dataset_id), ceiling(seq_along(unique(obis_list$dataset_id))/200))

pb <- progress::progress_bar$new(total = length(obis_ids))
obis_ds_context <- lapply(obis_ids, function(di) {
    dc <- robis::dataset(datasetid = di)
    dc |> 
        select(dataset_id = id, title, citation, description = abstract,
               nodes) |> 
        group_by(dataset_id) |> 
        mutate(nodes = paste(
            unlist(
                lapply(nodes, function(x){
                    paste(x[,"id"], x[,"name"], sep = ";")
                })
            ), collapse = "|"
        )) 
})

final_list <- bind_rows(obis_list, gbif_list)

gbif_ds_context <- bind_rows(gbif_ds_context)
obis_ds_context <- bind_rows(obis_ds_context)

gbif_ds_context$source <- "gbif"
obis_ds_context$source <- "obis"

final_context <- bind_rows(obis_ds_context, gbif_ds_context)

write_parquet(final_list, "data/reg_datasets_species.parquet")
write_parquet(final_context, "data/reg_datasets_context.parquet")

# Create and save simplified format
format_nodes <- function(node) {
    node <- strsplit(node, split = "\\|")[[1]]
    node <- unique(
        unlist(
            lapply(strsplit(node, split = ";"), function(x) {
                if (length(x) > 1) {
                    x[[2]]
                } else {
                    NULL
                }
            }),
            use.names = F
        )
    )
    paste(node, collapse = "|")
}
final_context_simp <- final_context |> 
    mutate(url = ifelse(
        source == "obis",
        paste0("https://obis.org/dataset/", dataset_id),
        paste0("https://www.gbif.org/dataset/", dataset_id)
    )) |> 
    group_by(title) |> 
    mutate(source = ifelse(length(source) > 1, paste(unique(source), collapse = "|"), source)) |> 
    distinct(title, .keep_all = TRUE) |> 
    select(Title = title, Source = source, Nodes = nodes, URL = url) |> 
    mutate(Source = toupper(Source)) |> 
    mutate(Nodes = format_nodes(Nodes))

jsonlite::write_json(list(
    project = "MPA Europe",
    creator = "OBIS",
    version = version,
    date = cur_date,
    datasets = final_context_simp
), "datasets_citation.json", pretty = TRUE)


# md_tab <- knitr::kable(final_context_simp)

# writeLines(
#     c(
#         "# Datasets used in MPA Europe SDMs",
#         "This work was conducted using data from OBIS and GBIF, encompassing thousands of datasets. 
#         We acknowledge the essential contributions of data providers who have made this information 
#         openly available through these platforms, as well as the extensive efforts of nodes in 
#         ensuring the data flow into the systems adheres to interoperable standards. 
#         Below is a list of all datasets utilized in this work. By following the provided links, 
#         you can find detailed information on how to cite each dataset. \n\n Note that there might be some duplicated
#         datasets, in case different names were used in both databases. 
#         Such duplicated data does not affect our analysis, since
#         this is filtered in the data-processing step. \n\n If a dataset is tagged in `Source` as `OBIS|GBIF`, that means that the dataset was identified both in OBIS and GBIF.
#         In that case, we provide the link for the OBIS version.", "",
#         md_tab
#     ),
#     con = "datasets.md"
# )

# END