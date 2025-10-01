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

get_done_species <- function(results_folder) {
    species <- list.files(results_folder)
    species <- gsub("taxonid=", "", species)
    return(species)
}

# Create a function to check good model -------
check_model <- function(results_folder, sp) {
    l <- jsonlite::read_json(glue(
        "{results_folder}/taxonid={sp}/model=mpaeu/taxonid={sp}_model=mpaeu_what=log.json"
    ))
    goodm <- unlist(l$model_good, use.names = F)
    if (is.numeric(goodm)) {
        goodm <- "esm"
    }
    if ("ensemble" %in% goodm) {
        best_m <- "ensemble"
    } else {
        evals <- lapply(goodm, \(mm) {
            mmb <- ifelse(mm == "rf", "rf_classification_ds", mm)
            eval <- arrow::read_parquet(glue(
                "{results_folder}/taxonid={sp}/model=mpaeu/metrics/taxonid={sp}_model=mpaeu_method={mmb}_what=cvmetrics.parquet"
            ))
            eval <- eval |>
                select(AUC = auc, CBI = cbi, TSS = tss_maxsss)
            eval <- round(apply(eval, 2, mean, na.rm = T), 2)
            eval <- as.data.frame(t(eval))
            eval$model <- mm
            eval
        })
        evals <- bind_rows(evals)

        best_m <- evals$model[which.max(evals$CBI)]
        best_m <- ifelse(best_m == "rf", "rf_classification_ds", best_m)
    }
    return(best_m)
}

# Create a function to do the processing --------
proc_layers <- function(sp, results_folder, out_folder, global_mask, base_file = "same",
                        sel_threshold = "p10", type = "std",
                        model_acro = "mpaeu", max_dispersal_const = 100000,
                        check_exists = TRUE, verbose = TRUE, max_mem = NULL) {

    if (!is.null(max_mem)) terra::terraOptions(memfrac = max_mem)
    if (verbose) message("Processing ", sp)

    if (!type %in% c("std", "const")) stop("`type` should be one of 'std' or 'const'")

    best_m <- check_model(results_folder, sp)

    global_mask <- vect(global_mask)

    # Load mask
    masks <- rast(glue(
        "{results_folder}/taxonid={sp}/model={model_acro}/predictions/taxonid={sp}_model={model_acro}_what=mask_cog.tif"
    ))
    masks <- subset(masks, ifelse(type == "std", "fit_region", "fit_region_max_depth"))
    NAflag(masks) <- 0

    # Load threshold
    thresholds <- arrow::read_parquet(glue(
        "{results_folder}/taxonid={sp}/model={model_acro}/metrics/taxonid={sp}_model={model_acro}_what=thresholds.parquet"
    ))

    th <- thresholds[grepl(strtrim(best_m, 2), thresholds$model), ]
    if (sel_threshold == "mss") {
        th <- th |> pull(max_spec_sens) * 100
    } else if (sel_threshold == "p10") {
        th <- th |> pull(p10) * 100
    } else {
        stop("Threshold not available")
    }

    all_preds <- list.files(glue("{results_folder}/taxonid={sp}/model={model_acro}/predictions"), full.names = T)
    all_preds <- all_preds[grepl(best_m, all_preds)]

    model_preds <- all_preds[!grepl("what=bootcv", all_preds)]

    if (base_file == "same") {
        base <- model_preds[1] |>
            rast()
    } else {
        base <- base_file |>
            rast()
    }

    base <- base |>
        terra::classify(matrix(data = c(-Inf, Inf, 0), nrow = 1), right = FALSE) |>
        crop(global_mask) |>
        mask(global_mask)

    if (type == "const") {
        # Load points and apply buffer
        fit_pts <- arrow::read_parquet(file.path(
            results_folder,
            glue::glue("taxonid={sp}/model={model_acro}/taxonid={sp}_model={model_acro}_what=fitocc.parquet")
        ))
        fit_pts <- terra::vect(as.data.frame(fit_pts), geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
        fit_pts_buffer <- terra::buffer(fit_pts, max_dispersal_const)
    }

    if (verbose) {
        pb <- progress::progress_bar$new(total = length(model_preds))
    }

    for (lp in model_preds) {

        if (verbose) pb$tick()

        basef <- basename(lp)
        basef <- gsub("_cog.tif", glue("_th={sel_threshold}_type={type}.tif"), basef)
        if (check_exists) {
            if (file.exists(file.path(out_folder, basef))) next
        }

        lyr_pred <- rast(lp)[[1]]

        lyr_pred <- terra::classify(lyr_pred, matrix(data = c(-Inf, th, 0), nrow = 1), right = FALSE)

        lyr_pred <- lyr_pred |>
            mask(masks) |>
            mask(global_mask) |>
            crop(global_mask) |>
            sum(base, na.rm = TRUE)

        if (all(is.na(minmax(lyr_pred)))) {
            return(paste0(sp, "_empty-on-starea"))
        } else if (minmax(lyr_pred)["min", ] < 0) {
            lyr_pred <- terra::classify(lyr_pred, matrix(data = c(-Inf, 0, 0), nrow = 1), right = FALSE)
        }
        if (max(minmax(lyr_pred)[,1]) == 0) {
            return(paste0(sp, "_all-zero-on-starea"))
        }

        if (type == "const") {
            pred_masked <- terra::mask(lyr_pred, lyr_pred != 0, maskvalues = 0)

            # Label patches and get polygon
            pred_patches <- terra::patches(pred_masked, directions = 8)
            patches_pols <- terra::as.polygons(pred_patches, dissolve = TRUE)

            # Find intersections
            touched_pt <- terra::intersect(patches_pols, fit_pts_buffer)
            touched_pt_ids <- unique(terra::values(touched_pt)[[1]])

            # Mask original
            touched_pt_f <- terra::classify(pred_patches, rcl = cbind(setdiff(unique(terra::values(pred_patches)), touched_pt_ids), NA))
            touched_pt_f <- terra::classify(touched_pt_f, matrix(c(-Inf, Inf, 1), ncol = 3))
            lyr_pred <- pred_masked |>
                mask(touched_pt_f) |>
                sum(base, na.rm = TRUE)

            if (minmax(lyr_pred)["min", ] < 0) {
                lyr_pred <- terra::classify(lyr_pred, matrix(data = c(-Inf, 0, 0), nrow = 1), right = FALSE)
            }
        }

        writeRaster(
            lyr_pred,
            file.path(out_folder, basef),
            datatype = "INT1U",
            overwrite = TRUE
        )
    }

    return(paste0(sp, "_done"))
}

prepare_layers <- function(
    thresholds = c("p10", "mss"),
    type = c("std", "const"),
    results_folder,
    out_folder,
    parallel = FALSE,
    n_cores = 80,
    max_mem = NULL, # use TRUE to use standard method (0.9/n_cores)
    global_mask = "data/shapefiles/mpa_europe_starea_v3.gpkg",
    base_file = "data/env/current/thetao_baseline_depthsurf_mean.tif",
    species = NULL, # pass species here to speed up, otherwise listed from folder
    verbose = FALSE
) {
    require(glue)
    require(dplyr)
    require(terra)

    fs::dir_create(out_folder)
    
    if (is.logical(max_mem) && isTRUE(max_mem)) max_mem <- (0.9 / n_cores)

    if (is.null(species)) species <- get_done_species()

    # Apply processing -----
    post_grid <- expand.grid(
        sel_threshold = thresholds,
        type = type
    )

    if (parallel) {
        require(furrr)
        plan(multisession, workers = n_cores)

        for (g in seq_len(nrow(post_grid))) {
            message("Processing ", post_grid$sel_threshold[g], " | ", post_grid$type[g], " (", g, ")\n")
            proc_result <- future_map(species, proc_layers,
                results_folder, out_folder, global_mask,
                sel_threshold = post_grid$sel_threshold[g],
                type = post_grid$type[g],
                model_acro = "mpaeu", base_file = base_file,
                check_exists = TRUE, verbose = FALSE,
                max_mem = max_mem, .progress = verbose
            )
            gc()
        }

        plan(sequential)
    } else {
        for (g in seq_len(nrow(post_grid))) {
            message("Processing ", post_grid$sel_threshold[g], " | ", post_grid$type[g], " (", g, ")\n")
            proc_result <- lapply(species, proc_layers,
                results_folder, out_folder, global_mask,
                sel_threshold = post_grid$sel_threshold[g],
                type = post_grid$type[g],
                model_acro = "mpaeu", base_file = base_file,
                check_exists = TRUE, verbose = verbose
            )
            gc()
        }
    }

    if (verbose) message("Processing concluded.")
    return(invisible(NULL))
}

list_processed <- function(out_folder, write_csv = TRUE, csv_folder = NULL) {
  if (is.null(csv_folder)) csv_folder <- getwd()

  done_files <- list.files(out_folder)
  done_files_m <- gsub("\\.tif", "", done_files)
  done_files_m <- gsub("rf_classification_ds", "rf", done_files_m)
  done_files_m <- gsub("current", "current_dec0", done_files_m)

  split_parts <- data.table::tstrsplit(done_files_m, "_")
  split_parts <- unlist(split_parts, use.names = F)
  split_parts <- gsub(".*=", "", split_parts)

  general_files <- as.data.frame(matrix(split_parts, nrow = length(done_files), byrow = FALSE))
  colnames(general_files) <- c("taxonid", "model", "method", "scenario", "period", "threshold", "type")
  general_files$file <- done_files

  general_files$period[general_files$period == "dec0"] <- NA
  general_files$method[general_files$method == "rf"] <- "rf_classification_ds"

  if (write_csv) {
    write.csv(general_files, file.path(csv_folder, "processed_files.csv"), row.names = FALSE)
    return(invisible(NULL))
  } else {
    return(general_files)
  }
}

zip_processed <- function(out_folder, by_type = TRUE, save_folder = NULL) {
  done_files <- list.files(out_folder)
  if (is.null(save_folder)) save_folder <- getwd()

  zip_files <- c()

  if (by_type) {
    p10_std <- done_files[grepl("th=p10_type=std", done_files)]
    mss_std <- done_files[grepl("th=mss_type=std", done_files)]
    p10_cons <- done_files[grepl("th=p10_type=cons", done_files)]
    mss_cons <- done_files[grepl("th=mss_type=cons", done_files)]

    if (length(p10_std) > 0) {
      message("Zipping p10 std")
      zip::zip(
        file.path(save_folder, "processed_p10_std.zip"),
        file.path(out_folder, p10_std),
        mode = "cherry-pick"
      )
      zip_files <- c(zip_files, file.path(save_folder, "processed_p10_std.zip"))
    }
    if (length(mss_std) > 0) {
      message("Zipping mss std")
      zip::zip(
        file.path(save_folder, "processed_mss_std.zip"),
        file.path(out_folder, mss_std),
        mode = "cherry-pick"
      )
      zip_files <- c(zip_files, file.path(save_folder, "processed_mss_std.zip"))
    }
    if (length(p10_cons) > 0) {
      message("Zipping p10 cons")
      zip::zip(
        file.path(save_folder, "processed_p10_cons.zip"),
        file.path(out_folder, p10_cons),
        mode = "cherry-pick"
      )
      zip_files <- c(zip_files, file.path(save_folder, "processed_mss_std.zip"))
    }
    if (length(mss_cons) > 0) {
      message("Zipping mss cons")
      zip::zip(
        file.path(save_folder, "processed_mss_cons.zip"),
        file.path(out_folder, mss_cons),
        mode = "cherry-pick"
      )
      zip_files <- c(zip_files, file.path(save_folder, "processed_mss_cons.zip"))
    }
  } else {
    message("Zipping all - this may take some time...")
    zip::zip(
      file.path(save_folder, "processed_all.zip"),
      file.path(out_folder, done_files),
      mode = "cherry-pick"
    )
    zip_files <- c(zip_files, file.path(save_folder, "processed_all.zip"))
  }
  cli::cli_alert_success("Files saved at {.file {zip_files}}")
  return(invisible(NULL))
}

# SESAM richness correction
# Check https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.13608
Rcpp::cppFunction('
NumericVector sesam_fast(NumericVector cell) {
  // Converted from R code with help of ChatGPT
  // first element is the richness
  int rich = (int)round(cell[0]);
  // The rest is the probabilities
  int nsp = cell.size() - 1;
  NumericVector probs(nsp);
  for (int i = 0; i < nsp; i++) {
    probs[i] = cell[i + 1];
  }
  // If any is NA, need to return NA
  for (int i = 0; i < nsp; i++) {
    if (NumericVector::is_na(probs[i]) || NumericVector::is_na(cell[0])) {
      NumericVector out(nsp, NA_REAL);
      return out;
    }
  }
  // If not proceed
  if (rich > 0) {
    // vector to hold indices
    std::vector<int> idx(nsp);
    for (int i = 0; i < nsp; i++) idx[i] = i;

    // sort indices by probs (descending)
    std::sort(idx.begin(), idx.end(),
      [&](int a, int b) { return probs[a] > probs[b]; });
    // convert to 0/1
    NumericVector out(nsp);
    for (int i = 0; i < rich && i < nsp; i++) {
      out[idx[i]] = 1.0;
    }
    // others remain 0
    return out;
  } else {
    // No richness, returns 0
    NumericVector out(nsp);
    std::fill(out.begin(), out.end(), 0.0);
    return out;
  }
}
')

sesam_prr <- function(probabilities, richness) {
    # Modified version of ecospat::ecospat.SESAM.prr to work with terra raster
    # terra::app(c(richness, probabilities), \(cell) {
    #     nc <- names(cell)[-1]
    #     cell <- unlist(cell, use.names = F)
    #     rich <- as.integer(round(cell[1]))
    #     probs <- cell[2:length(cell)]
    #     if (rich > 0) {
    #         com <- order(probs, decreasing = TRUE)
    #         pres <- com[1:rich]
    #         probs[pres] <- 1
    #         probs[-pres] <- 0
    #     } else {
    #         probs[] <- 0
    #     }
    #     names(probs) <- nc
    #     return(probs)
    # })
    terra::app(c(richness, probabilities), sesam_fast)
}

# Validation of this function
# library(ecospat)
# proba <- ecospat.testData[,73:92]
# sr <- as.data.frame(rowSums(proba))
# ppr<-ecospat.SESAM.prr(proba, sr)
# head(ppr)
# # Test ours
# library(terra)
# prob_r <- rast(nlyr = 20, ncol = 10, nrow = 30)
# for (i in 1:20) {
#     prob_r[[i]][] <- proba[,i]
# }
# ric_r <- sum(prob_r)
# ppr_v2 <- sesam_prr(prob_r, ric_r)
# ppr <- rowSums(ppr)
# ppr_v2 <- sum(ppr_v2)
# all.equal(as.vector(ppr), as.vector(ppr_v2[]))


# Implements a fast C++ code of the LCBD as done in adespatial::beta.div()
Rcpp::cppFunction('
NumericVector lcbd_cpp(const NumericMatrix& Y) {
  // Generated using LLM tool: ChatGPT

  int n = Y.nrow();
  int m = Y.ncol();
  
  // Hellinger transform
  NumericMatrix Yh(n, m);
  for (int i = 0; i < n; i++) {
    double rowsum = 0.0;
    for (int j = 0; j < m; j++) {
      rowsum += Y(i, j);
    }
    if (rowsum > 0) {
      for (int j = 0; j < m; j++) {
        Yh(i, j) = sqrt(Y(i, j) / rowsum);
      }
    }
  }
  
  // Compute column means (centroid)
  NumericVector centroid(m);
  for (int j = 0; j < m; j++) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
      sum += Yh(i, j);
    }
    centroid[j] = sum / n;
  }
  
  // Compute squared distances to centroid
  NumericVector d2(n);
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    double dist2 = 0.0;
    for (int j = 0; j < m; j++) {
      double diff = Yh(i, j) - centroid[j];
      dist2 += diff * diff;
    }
    d2[i] = dist2;
    total += dist2;
  }
  
  // LCBD values
  NumericVector LCBD(n);
  for (int i = 0; i < n; i++) {
    LCBD[i] = d2[i] / total;
  }
  
  return LCBD;
}
')

# Validate code
# library(adespatial)
# library(vegan)
# data(mite)
# tictoc::tic()
# res <- beta.div(mite, adj = F)
# tictoc::toc()
# tictoc::tic()
# res_cpp <- lcbd_cpp(as.matrix(mite))
# tictoc::toc()
# all.equal(unname(res$LCBD), res_cpp)







# Old code, for back compatibility - to be removed soon

raw_processing <- function(gr, sp_list, base_raw, study_area, output_folder, acro, eez_cell = NULL, protseas_cell = NULL) {
  
  grsp_raw <- sp_list$taxonID[sp_list$group == gr]
  #groups_sp[[gr]][["raw"]] <- grsp_raw
  grsp_raw <- paste0("data/species/key=", grsp_raw, ".parquet")
  grsp_raw <- grsp_raw[file.exists(grsp_raw)]
  
  cli::cli_alert_info("{length(grsp_raw)} species available for raw based join for group {gr}. Joining...")
  base_raw <- terra::rast(base_raw)
  study_area <- terra::vect(study_area)

  base_raw <- terra::classify(base_raw, cbind(-100, 100, 0))

  gr_raw <- base_raw
  gr_raw_agg <- terra::aggregate(base_raw, fact = 10, na.rm = T)

  base_raw <- terra::mask(base_raw, study_area)
  
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

  gr_raw <- terra::mask(gr_raw, base_raw)
  base_raw <- terra::aggregate(base_raw, fact = 10, na.rm = T)
  gr_raw_agg <- terra::mask(gr_raw_agg, base_raw)
  
  outraw <- glue::glue(
    "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_type=original.tif"
  )
  outrawagg <- glue::glue(
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
  
  # # EEZ and ProtectedSeas
  # outeez <- glue(
  #   "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_area=eez.txt"
  # )
  
  # outprot <- glue(
  #   "metric=richness_model={acro}_method=raw_scen=current_group={gsub('/', '-', tolower(gr))}_area=mpa.txt"
  # )
  
  # pt$cell <- NA
  # pt$cell <- cellFromXY(gr_raw, as.data.frame(pt[, c("decimalLongitude", "decimalLatitude")]))
  
  # pt_eez <- left_join(eez_cell, pt, relationship = "many-to-many")
  
  # pt_eez_ed <- pt_eez %>%
  #   filter(!is.na(taxonID)) %>%
  #   group_by(MRGID) %>%
  #   distinct(taxonID) %>%
  #   summarise(total = n())
  
  # write.table(pt_eez_ed, file.path(output_folder, outeez), row.names = F)
  
  # pt_prot <- left_join(protseas_cell, pt, relationship = "many-to-many")
  
  # pt_prot_ed <- pt_prot %>%
  #   filter(!is.na(taxonID)) %>%
  #   group_by(SITE_ID) %>%
  #   distinct(taxonID) %>%
  #   summarise(total = n())
  
  # write.table(pt_prot_ed, file.path(output_folder, outprot), row.names = F)
  
  # return(invisible(NULL))
  return(file.path(output_folder, outraw))
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
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_what=mask_cog.tif"
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
      "taxonid={ids}/model=mpaeu/predictions/taxonid={ids}_model=mpaeu_what=mask_cog.tif"
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

  if (length(model) > 1) {
    model <- "combined"
  }

  if (save_parquet) {
    names(binary) <- paste0("taxonid=", ids)
    binary_df <- as.data.frame(binary, xy = T)
    #binary_df <- tidyr::pivot_longer(binary_df, 3:ncol(binary_df), names_to = "taxonid", values_to = "value")
    arrow::write_parquet(
      binary_df,
      tolower(file.path(outfolder, paste0("richness_", group, "_", model, "_", scenario, "_binary_part", index, ".parquet")))
    )
    rm(binary_df)
  }

  # Sum results
  classified <- classified/100
  classified <- sum(classified, na.rm = TRUE)

  binary <- sum(binary, na.rm = TRUE)

  # Save results
  outf <- tolower(file.path(outfolder, paste0("richness_", group, "_", model, "_", scenario, "_binary_part", index, ".tif")))

  terra::writeRaster(binary, outf, overwrite = TRUE)
  terra::writeRaster(classified, gsub("binary", "cont", outf), overwrite = TRUE)

  terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE)

  return(outf)
}
