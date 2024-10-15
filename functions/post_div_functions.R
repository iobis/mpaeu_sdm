############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
####################### Diversity maps - functions #############################

get_sp_list <- function() {
    sp_list <- obissdm::recent_file("data", "all_splist")
    sp_list <- read.csv(sp_list)

    sp_list <- sp_list %>%
        mutate(group = case_when(
            kingdom == "Chromista" ~ "Chromista",
            kingdom == "Plantae" ~ "Plantae",
            phylum == "Annelida" ~ "Annelida",
            phylum == "Arthropoda" ~ "Arthropoda",
            phylum == "Cnidaria" ~ "Cnidaria",
            phylum == "Echinodermata" ~ "Echinodermata",
            phylum == "Mollusca" ~ "Mollusca",
            phylum == "Nematoda" ~ "Nematoda",
            phylum == "Chordata" ~ "Chordata",
            .default = "Others"
        )) %>%
        mutate(group = case_when(
            class == "Aves" ~ "Aves",
            class == "Mammalia" ~ "Mammalia",
            class == "Myxini" | class == "Petromyzonti" ~ "Myxini/Petromyzonti",
            class == "Teleostei" | class == "Elasmobranchii" |
                class == "Holocephali" | class == "Chondrostei" ~ "Fishes",
            .default = group
        ))

    return(sp_list)
}

raw_processing <- function(sp_list, gr, base_raw, output_folder) {

    grsp_raw <- sp_list$taxonID[sp_list$group == gr]
    groups_sp[[gr]][["raw"]] <- grsp_raw
    grsp_raw <- paste0("data/species/key=", grsp_raw, ".parquet")
    grsp_raw <- grsp_raw[file.exists(grsp_raw)]

    cli::cli_alert_info("{length(grsp_raw)} species available for raw based join. Joining...")
    gr_raw <- base_raw
    gr_raw_agg <- terra::aggregate(base_raw, fact = 10)

    pt <- arrow::open_dataset(grsp_raw)
    pt <- pt %>%
        filter(data_type == "fit_points") %>%
        select(decimalLongitude, decimalLatitude, taxonID) %>%
        collect()

    pt$cell <- cellFromXY(gr_raw, as.data.frame(pt[, c("decimalLongitude", "decimalLatitude")]))
    pt_ed <- pt %>%
        group_by(cell) %>%
        distinct(taxonID) %>%
        summarise(total = n())

    gr_raw[pt_ed$cell] <- pt_ed$total

    pt$cell <- NA
    pt$cell <- cellFromXY(gr_raw_agg, as.data.frame(pt[, c("decimalLongitude", "decimalLatitude")]))
    pt_ed <- pt %>%
        group_by(cell) %>%
        distinct(taxonID) %>%
        summarise(total = n())

    gr_raw_agg[pt_ed$cell] <- pt_ed$total

    gr_raw <- as.int(gr_raw)
    gr_raw_agg <- as.int(gr_raw_agg)

    outraw <- glue(
        "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_type=original.tif"
    )
    outrawagg <- glue(
        "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_type=aggregated.tif"
    )

    if (terra::minmax(gr_raw_agg)[2, 1] <= 255) {
        format_out <- "INT1U"
    } else if (terra::minmax(gr_raw_agg)[2, 1] <= 65535) {
        format_out <- "INT2U"
    } else {
        format_out <- "INT4U"
    }

    terra::writeRaster(gr_raw, file.path(output_folder, outraw), overwrite = T, datatype = format_out)
    obissdm::cogeo_optim(file.path(output_folder, outraw))

    terra::writeRaster(gr_raw_agg, file.path(output_folder, outrawagg), overwrite = T, datatype = format_out)
    obissdm::cogeo_optim(file.path(output_folder, outrawagg))

    # EEZ and WDPA
    pt$cell <- NA
    pt$cell <- cellFromXY(gr_raw_agg, as.data.frame(pt[, c("decimalLongitude", "decimalLatitude")]))
    pt_ed <- pt %>%
        group_by(cell) %>%
        distinct(taxonID) %>%
        summarise(total = n())

    gr_raw_agg[pt_ed$cell] <- pt_ed$total

    return(invisible(NULL))
}


extract_shapes <- function(eez, wdpa, layers) {
    eez_extracted
    wdpa_extracted
}