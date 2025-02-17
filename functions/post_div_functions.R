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

raw_processing <- function(sp_list, gr, base_raw, output_folder, eez_cell, protseas_cell) {

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

    # EEZ and ProtectedSeas
    outeez <- glue(
      "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_area=eez.txt"
    )
    
    outprot <- glue(
      "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_area=mpa.txt"
    )
    
    pt$cell <- NA
    pt$cell <- cellFromXY(gr_raw, as.data.frame(pt[, c("decimalLongitude", "decimalLatitude")]))
    
    pt_eez <- left_join(eez_cell, pt, relationship = "many-to-many")
    
    pt_eez_ed <- pt_eez %>%
      filter(!is.na(taxonID)) %>%
      group_by(MRGID) %>%
      distinct(taxonID) %>%
      summarise(total = n())
    
    write.table(pt_eez_ed, file.path(output_folder, outeez), row.names = F)
    
    pt_prot <- left_join(protseas_cell, pt, relationship = "many-to-many")
    
    pt_prot_ed <- pt_prot %>%
      filter(!is.na(taxonID)) %>%
      group_by(SITE_ID) %>%
      distinct(taxonID) %>%
      summarise(total = n())
    
    write.table(pt_prot_ed, file.path(output_folder, outprot), row.names = F)

    return(invisible(NULL))
}


extract_shapes <- function(eez_cell, protseas_cell, layers) {
  
  # results_eez <- lapply(1:nlyr(layers), function(x) NULL)
  # results_prot <- results_eez
  
  total <- terra::nlyr(layers)
  
  layer_sources <- sources(layers)
  
  cl <- parallel::makeCluster(30)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, "layer_sources", envir=environment())
  results_general <- foreach(i = seq_len(total)) %dopar% {
    
    require("dplyr")
    
    lay <- terra::rast(layer_sources[i])
    
    eez_cell <- readRDS("data/shapefiles/eez_cell.rds")
    protseas_cell <- readRDS("data/shapefiles/protseas_cell.rds")
    
    l_eez <- lay[eez_cell$cell]
    l_prot <- lay[protseas_cell$cell]
    
    colnames(l_eez) <- colnames(l_prot) <- "presence"
    
    l_eez <- l_eez %>%
      mutate(MRGID = eez_cell$MRGID) %>%
      filter(!is.na(presence)) %>%
      group_by(MRGID) %>%
      summarise(presence = ifelse(any(presence == 1), i, 0)) %>%
      filter(presence != 0)
    
    l_prot <- l_prot %>%
      mutate(SITE_ID = protseas_cell$SITE_ID) %>%
      filter(!is.na(presence)) %>%
      group_by(SITE_ID) %>%
      summarise(presence = ifelse(any(presence == 1), i, 0)) %>%
      filter(presence != 0)
    
    list(eez = l_eez, mpa = l_prot)
  }
  parallel::stopCluster(cl)
  
    # for (i in 1:total) {
    #   cli::cli_alert_info("Processing layer {i} out of {total}")
    #   l_eez <- layers[[i]][eez_cell$cell]
    #   l_prot <- layers[[i]][protseas_cell$cell]
    #   
    #   colnames(l_eez) <- colnames(l_prot) <- "presence"
    #   
    #   l_eez <- l_eez %>%
    #     mutate(MRGID = eez_cell$MRGID) %>%
    #     filter(!is.na(presence)) %>%
    #     group_by(MRGID) %>%
    #     summarise(presence = ifelse(any(presence == 1), i, 0)) %>%
    #     filter(presence != 0)
    #   
    #   l_prot <- l_prot %>%
    #     mutate(SITE_ID = protseas_cell$SITE_ID) %>%
    #     filter(!is.na(presence)) %>%
    #     group_by(SITE_ID) %>%
    #     summarise(presence = ifelse(any(presence == 1), i, 0)) %>%
    #     filter(presence != 0)
    #   
    #   results_eez[[i]] <- l_eez
    #   results_prot[[i]] <- l_prot
    # }
  
  results_eez <- lapply(results_general, function(x) x$eez)
  results_prot <- lapply(results_general, function(x) x$mpa)
  
  results_eez <- dplyr::bind_rows(results_eez)
  results_prot <- dplyr::bind_rows(results_prot)
  
  results_eez <- results_eez %>%
    group_by(MRGID) %>%
    distinct(presence) %>%
    summarise(total = n())
  
  results_prot <- results_prot %>%
    group_by(SITE_ID) %>%
    distinct(presence) %>%
    summarise(total = n())
  
  return(list(eez = results_eez, mpa = results_prot))
  
}



# Convert polygons to cell
get_eez_cell <- function(eez, resolution = 0.05) {
  if (!file.exists("data/shapefiles/eez_cell.rds")) {
    
    eez_cell <- lapply(1:nrow(eez), function(x) NULL)
    r <- terra::rast(res = resolution)
    r[] <- 1
    
    for (i in 1:nrow(eez)) {
      cat(i, "\n")
      cell_r <- terra::extract(r, eez[i,], cell = T)
      eez_cell[[i]] <- data.frame(cell = cell_r$cell)
    }
    
    names(eez_cell) <- eez$MRGID
    
    eez_cell_bind <- dplyr::bind_rows(eez_cell, .id = "MRGID")
    
    saveRDS(eez_cell_bind, "data/shapefiles/eez_cell.rds")
  } else {
    eez_cell_bind <- readRDS("data/shapefiles/eez_cell.rds")
  }
  return(eez_cell_bind)
}

get_protseas_cell <- function(protseas, resolution = 0.05) {
  if (!file.exists("data/shapefiles/protseas_cell.rds")) {
    #protseas <- protseas[protseas$lfp > 2,"SITE_ID"]
    
    protseas_cell <- lapply(1:nrow(protseas), function(x) NULL)
    r <- terra::rast(res = resolution)
    r[] <- 1
    
    for (i in 1:nrow(protseas)) {
      cat(i, "\n")
      cell_r <- terra::extract(r, protseas[i,], cell = T)
      protseas_cell[[i]] <- data.frame(cell = cell_r$cell)
    }
    
    names(protseas_cell) <- protseas$SITE_ID
    
    protseas_cell_bind <- dplyr::bind_rows(protseas_cell, .id = "SITE_ID")
    
    saveRDS(protseas_cell_bind, "data/shapefiles/protseas_cell.rds")
  } else {
    protseas_cell_bind <- readRDS("data/shapefiles/protseas_cell.rds")
  }
  return(protseas_cell_bind)
}

# # Convert polygons to h3
# get_eez_h3 <- function(eez) {
#   if (!file.exists("data/shapefiles/eez_h3.rds")) {
#     eez_sf <- sf::st_as_sf(eez[,1])
#     
#     eez_sf_h3 <- lapply(1:nrow(eez_sf), function(x) NULL)
#     
#     for (i in 1:nrow(protseas_sf)) {
#       cat(i, "\n")
#       eez_sf_h3[[i]] <- data.frame(h3 = unlist(h3jsr::polygon_to_cells(eez_sf[i,], res = 7), use.names = F))
#     }
#     
#     names(eez_sf_h3) <- eez_sf$MRGID
#     
#     eez_sf_bind <- dplyr::bind_rows(eez_sf_h3, .id = "MRGID")
#     
#     saveRDS(eez_sf_bind, "data/shapefiles/eez_h3.rds")
#   } else {
#     eez_sf_bind <- readRDS("data/shapefiles/eez_h3.rds")
#   }
#   return(eez_sf_bind)
# }
# 
# get_protseas_h3 <- function(protseas) {
#   if (!file.exists("data/shapefiles/protseas_h3.rds")) {
#     protseas <- protseas[protseas$lfp > 2,]
#     protseas_sf <- sf::st_as_sf(protseas[,1])
#     
#     protseas_sf_h3 <- lapply(1:nrow(protseas_sf), function(x) NULL)
#     
#     for (i in 1:nrow(protseas_sf)) {
#       cat(i, "\n")
#       protseas_sf_h3[[i]] <- data.frame(h3 = unlist(h3jsr::polygon_to_cells(protseas_sf[i,], res = 7), use.names = F))
#     }
#     
#     names(protseas_sf_h3) <- protseas_sf$SITE_ID
#     
#     protseas_sf_bind <- dplyr::bind_rows(protseas_sf_h3, .id = "MRGID")
#     
#     saveRDS(protseas_sf_bind, "data/shapefiles/protseas_h3.rds")
#   } else {
#     protseas_sf_bind <- readRDS("data/shapefiles/protseas_h3.rds")
#   }
#   return(protseas_sf_bind)
# }


### Main richness function 
# You will need this:
# xr <- import("xarray")
# da <- import("dask")
# ri <- import("rioxarray")
# np <- import("numpy")
# dd <- import("dask.distributed")
# client <- dd$Client()
# browseURL(client$dashboard_link)


# Create function for processing
proc_maps_py <- function(index, batches, model, threshold_table, outfolder, group, scenario, results_folder, gen_cont = FALSE) {
  ids <- batches[[index]]
  
  if (!is.list(model)) {
    sel_threshold <- threshold_table[match(ids, threshold_table$taxonID, nomatch = 0), ]
    if (!all.equal(sel_threshold$taxonID, ids)) stop("Problem with thresholds IDs")
    
    sel_threshold <- pull(sel_threshold, model)
    sel_threshold <- sel_threshold[!is.na(sel_threshold)]
    if (length(sel_threshold) != length(ids)) stop("Problem with thresholds table")
  } else {
    model <- model[[index]]
    if (length(ids) != length(model)) stop("Problem in valid models vector")
    
    sel_threshold <- threshold_table[match(ids, threshold_table$taxonID, nomatch = 0), ]
    if (!all.equal(sel_threshold$taxonID, ids)) stop("Problem with thresholds IDs")
    sel_threshold <- unlist(lapply(1:nrow(sel_threshold), function(x) {
      pull(sel_threshold, model[x])[x]
    }))
  }
  
  modelf <- ifelse(model == "rf",
                   "rf_classification_ds", model
  )
  
  raster_paths <- file.path(
    results_folder,
    glue::glue(
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_method={modelf}_scen={scenario}_cog.tif"
    )
  )
  
  mask_paths <- file.path(
    results_folder,
    glue::glue(
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_mask_cog.tif"
    )
  )
  
  # Open xarray datasets
  rasters <- xr$open_mfdataset(
    raster_paths,
    engine = "rasterio", concat_dim = "band", combine = "nested", chunks = "auto"
  )
  
  masks <- xr$open_mfdataset(
    mask_paths,
    engine = "rasterio", concat_dim = "band", combine = "nested", chunks = "auto"
  )
  
  # Masks band 1 = native ecoregions
  masks_band1 <- masks$sel(band = 1)
  
  rasters_band1 <- rasters$sel(band = 1)
  
  masked_rasters <- rasters_band1 * masks_band1
  
  
  #### Produce binary version
  thresholds_da <- xr$DataArray(as.integer(as.vector(sel_threshold)),
                                dims = list("band"), coords = list("band" = masked_rasters$band)
  )
  
  # Apply thresholds: any value < threshold becomes 0, >= threshold becomes 1
  binary_rasters <- xr$where(
    masked_rasters$isnull(), # Check if values are NaN
    np$nan, # If NaN, keep as NaN
    xr$where(masked_rasters >= thresholds_da, 1, 0)
  )
  
  layer_sum <- binary_rasters$sum(dim = "band", skipna = FALSE)
  
  layer_sum_filled <- layer_sum$fillna(-1)
  
  layer_sum_int <- layer_sum_filled$astype("int32")
  
  
  #### Produce non-binary version
  if (gen_cont) {
    non_binary <- xr$where(
      masked_rasters$isnull(), # Check if values are NaN
      np$nan, # If NaN, keep as NaN
      xr$where(masked_rasters >= thresholds_da, masked_rasters, 0)
    )
    
    # Scale each raster layer to 0-1 by dividing by 100
    scaled_rasters <- non_binary / 100
    
    non_binary_layer_sum <- scaled_rasters$sum(dim = "band", skipna = FALSE)
    
    non_binary_layer_sum_filled <- non_binary_layer_sum$fillna(-1)
  }
  
  
  # Save results
  group <- gsub("\\/", "-", group)
  
  if (length(model) > 1) {
    model <- "combined"
  }
  
  outf <- file.path(outfolder, paste0("richness_", group, "_", model, "_", scenario, "_binary_part", index, ".tif"))
  layer_sum_int$rio$to_raster(outf)
  
  if (gen_cont) {
    non_binary_layer_sum_filled$rio$to_raster(gsub("binary", "cont", outf))
  }
  
  return(outf)
}



#### New version
proc_maps <- function(
    index,
    batches, # Recommended a maximum of 200 per batch
    model,
    threshold_table,
    outfolder,
    group,
    scenario,
    results_folder,
    studyarea,
    mask_layer,
    save_parquet = TRUE,
    max_mem = 0.06) {

  terra::terraOptions(memfrac = max_mem)
  
  ids <- batches[[index]]

  if (!is.list(model)) {
    sel_threshold <- threshold_table[match(ids, threshold_table$taxonID, nomatch = 0), ]
    if (!all.equal(sel_threshold$taxonID, ids)) stop("Problem with thresholds IDs")

    sel_threshold <- pull(sel_threshold, model)
    sel_threshold <- sel_threshold[!is.na(sel_threshold)]
    if (length(sel_threshold) != length(ids)) stop("Problem with thresholds table")
  } else {
    model <- model[[index]]
    if (length(ids) != length(model)) stop("Problem in valid models vector")

    sel_threshold <- threshold_table[match(ids, threshold_table$taxonID, nomatch = 0), ]
    if (!all.equal(sel_threshold$taxonID, ids)) stop("Problem with thresholds IDs")
    sel_threshold <- unlist(lapply(1:nrow(sel_threshold), function(x) {
      pull(sel_threshold, model[x])[x]
    }))
  }

  # Load C++ function
  Rcpp::sourceCpp("functions/classify_raster.cpp")

  modelf <- ifelse(model == "rf",
    "rf_classification_ds", model
  )

  studyarea <- terra::vect(studyarea)

  raster_paths <- file.path(
    results_folder,
    glue::glue(
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_method={modelf}_scen={scenario}_cog.tif"
    )
  )

  mask_paths <- file.path(
    results_folder,
    glue::glue(
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_mask_cog.tif"
    )
  )

  # Open layers
  rasters <- lapply(raster_paths, function(x){
    rr <- terra::rast(x, lyrs = 1)
    if (ext(rr) != ext(-180, 180, -90, 90)) {
      rr <- extend(rr, ext(-180, 180, -90, 90))
    }
    rr
  })
  rasters <- terra::rast(rasters)
  
  masks <- lapply(mask_paths, function(x){
    rr <- terra::rast(x, lyrs = mask_layer)
    if (ext(rr) != ext(-180, 180, -90, 90)) {
      rr <- extend(rr, ext(-180, 180, -90, 90))
    }
    rr
  })
  masks <- terra::rast(masks)
  terra::NAflag(masks) <- 0
  masks <- terra::mask(masks, studyarea)

  rasters <- terra::crop(rasters, studyarea)
  masks <- terra::crop(masks, studyarea)

  dsr <- terra::sds(list(rasters, masks))

  # Produce masked version
  masked <- terra::app(dsr, prod)

  # Produce continuous version
  classified <- terra::app(masked,
                           applyThreshold,
                           thresholds = as.integer(as.vector(sel_threshold)))

  # Produce binary version
  binary <- terra::classify(classified, cbind(0, Inf, 1))
  binary <- terra::as.int(binary)

  group <- gsub("\\/", "-", group)

  if (save_parquet) {
    names(binary) <- paste0("taxonid=", ids)
    binary_df <- as.data.frame(binary, xy = T)
    #binary_df <- tidyr::pivot_longer(binary_df, 3:ncol(binary_df), names_to = "taxonid", values_to = "value")
    arrow::write_parquet(
      binary_df,
      file.path(outfolder, paste0("richness_", group, "_", model, "_", scenario, "_binary_part", index, ".parquet"))
    )
    rm(binary_df)
  }

  # Sum results
  classified <- classified/100
  classified <- sum(classified, na.rm = TRUE)

  binary <- sum(binary, na.rm = TRUE)

  # Save results
  if (length(model) > 1) {
    model <- "combined"
  }

  outf <- file.path(outfolder, paste0("richness_", group, "_", model, "_", scenario, "_binary_part", index, ".tif"))

  terra::writeRaster(binary, outf, overwrite = TRUE)
  terra::writeRaster(classified, gsub("binary", "cont", outf), overwrite = TRUE)

  terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE)

  return(outf)
}
