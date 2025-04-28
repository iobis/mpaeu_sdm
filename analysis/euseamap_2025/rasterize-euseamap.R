# Rasterize EU-SeaMap
library(sf)
sf_use_s2(FALSE)
library(exactextractr)
library(terra)
library(dplyr)
outf <- "proc-layers"
fs::dir_create(outf)

# Path to the unziped EUSeaMap file
euseamap <- st_read("~/Downloads/EUSeaMap_2023/EUSeaMap_2023.gdb")

un_codes <- unique(euseamap$EUNIScomb)
extent_eusea <- ext(euseamap)

features_list <- lapply(un_codes, function(x) NULL)

pb <- progress::progress_bar$new(total = length(un_codes))

for (i in seq_len(length(un_codes))) {

    pb$tick()

    sel_layer <- euseamap[euseamap$EUNIScomb == un_codes[i],]

    sel_layer <- suppressMessages(
        suppressWarnings(
            sel_layer %>% group_by(EUNIS2019C) %>%
                summarise(geometry = st_union(Shape)) %>%
                ungroup() 
        )
    )

    features_list[[i]] <- sel_layer %>% select(EUNIS2019C) %>% st_drop_geometry()
    features_list[[i]]$inside_id <- seq_len(nrow(features_list[[i]]))
    features_list[[i]]$EUNIScomb <- un_codes[i]
    features_list[[i]]$proc_layer <- i

    if (any(st_geometry_type(sel_layer) == "POLYGON")) {
        sel_layer <- st_cast(sel_layer, "MULTIPOLYGON")
    }

    frac <- coverage_fraction(raster::raster(res = 0.05), sel_layer, crop = FALSE)
    frac <- lapply(frac, terra::rast)

    frac_all <- terra::rast(frac)
    frac_all[frac_all == 0] <- NA

    frac_all_max <- which.max(frac_all)

    frac_crop <- crop(frac_all_max, extent_eusea)
    frac_crop <- as.int(frac_crop)
    names(frac_crop) <- un_codes[i]

    writeRaster(frac_crop, file.path(outf, paste0("proc_layer_EUComb", i, ".tif")), overwrite = T)

}

features_list <- do.call("rbind", features_list)

all_list <- euseamap %>% select(
    EUNIScomb, EUNIScombD, EUNIS2019C, EUNIS2019D
) %>% st_drop_geometry() %>% distinct()

features_list <- left_join(features_list, all_list)

write.csv(features_list, "features_list_EUSeaMap.csv", row.names = F)
