############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# June of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################## Get thermal range maps ##############################

# Load packages ----
library(terra)
library(arrow)
library(dplyr)
library(furrr)
source("functions/components_model_species.R")
source("functions/utils.R")

# Settings ----
model_acro <- "mpaeu"
output_format <- "multiband" # either "parquet" or "multiband"
results_folder <- "/data/scps/v5/results"
tg_species <- as.integer(extract_storr())

# Prepare adapted function to get thermal envelope
# This function is based on the `get_thermal_envelope` of the package `speedy`
# but use the data we used for SDMs instead
# https://github.com/iobis/rspeedy

# This function to overcome error in findInterval
alt_kde <- function(p, fhat) {
  if (any(p > 1) | any(p < 0)) 
    stop("p must be <= 1 and >= 0")
  cumul.prob <- ks::pkde(q = fhat$eval.points, fhat = fhat)
  ind <- findInterval(x = p, vec = round(cumul.prob, 10))
  quant <- rep(0, length(ind))
  for (j in 1:length(ind)) {
    i <- ind[j]
    if (i == 0) 
      quant[j] <- fhat$eval.points[1]
    else if (i >= length(fhat$eval.points)) 
      quant[j] <- fhat$eval.points[length(fhat$eval.points)]
    else {
      quant1 <- fhat$eval.points[i]
      quant2 <- fhat$eval.points[i + 1]
      prob1 <- cumul.prob[i]
      prob2 <- cumul.prob[i + 1]
      alpha <- (p[j] - prob2)/(prob1 - prob2)
      quant[j] <- quant1 * alpha + quant2 * (1 - alpha)
    }
  }
  return(quant)
}

get_envelope <- function (sp_pts, env) {
 
  # Use same standards as speedy
  #sf::sf_use_s2(FALSE)
  #points_per_m2 <- 1/3e+08
  max_points <- 10000
  kde_bandwidth <- 0.8
  #min_area <- 2e+09
  min_temperatures <- 3
  
  dist <- vect(sp_pts, geom = c("decimalLongitude", "decimalLatitude"),
               crs = "EPSG:4326")
  
  if (nrow(dist) == 0) {
    return(NULL)
  }
  
  if (nrow(dist) > max_points) {
    points <- dist[sample(seq_len(nrow(dist)), max_points),]
  } else {
    points <- dist
  }
  
  temperatures <- terra::extract(env, points)
  temperatures <- temperatures[,2]
  temperatures <- temperatures[!is.na(temperatures)]
  
  if (length(temperatures) < min_temperatures) {
    return(NULL)
  }
  
  kd <- ks::kde(temperatures, h = kde_bandwidth)
  percentiles <- try(ks::qkde(c(0.01, 0.99), kd), silent = T)
  if (inherits(percentiles, "try-error")) {
    percentiles <- alt_kde(c(0.01, 0.99), kd)
  }
  
  return(percentiles)
}

get_masked <- function(percentiles, layer) {

  masked <- terra::classify(layer, matrix(c(-Inf, percentiles[1], NA), ncol = 3), right = F)
  masked <- terra::classify(masked, matrix(c(percentiles[2], Inf, NA), ncol = 3))
  masked <- terra::classify(masked, matrix(c(-Inf, Inf, 1), ncol = 3))
  
  return(masked)
}

# Create a function to generate the thermal range maps
get_thermrange <- function(species, target_folder, model_acro, output_format,
                           skip_done = TRUE, mem_frac = NULL) {
  
  # For memory management
  if (!is.null(mem_frac)) terra::terraOptions(memfrac = mem_frac)

  # Load species data
  sp_data <- read_parquet(paste0("data/species/key=", species, ".parquet"))
  
  sp_sel_data <- sp_data %>%
    select(decimalLongitude, decimalLatitude, data_type, taxonID) %>%
    filter(data_type == "fit_points")
  
  if (output_format == "multiband") {
    outfile <- paste0(target_folder, "/taxonid=", species, "/model=", model_acro,
                    "/predictions/taxonid=", species, "_model=",
                    model_acro, "_what=thermenvelope.tif")
  } else if (output_format == "parquet") {
    outfile <- paste0(target_folder, "/taxonid=", species, "/model=", model_acro,
                    "/predictions/taxonid=", species, "_model=",
                    model_acro, "_what=thermenvelope.parquet")
  } else {
    cli::cli_abort("`output_format` not recognized.")
  }
  
  if (file.exists(outfile) & skip_done) {
    return("already_saved")
  }
  
  if (nrow(sp_sel_data) >= 10) {
    
    # Load ecological information
    eco_info <- arrow::open_csv_dataset("data/species_ecoinfo.csv") %>%
      filter(taxonID == species) %>%
      collect()
    
    if (nrow(eco_info) < 1) {
      eco_info <- obissdm::mp_get_ecoinfo(
        species_list = species,
        outfile = NULL,
        return_table = TRUE,
        show_progress = FALSE
      )
    }
    
    # Select ecological information for the species
    hab <- eco_info$mode_life
    hab_depth <- hab_to_depth(hab)
    
    # Load temperature data
    temp_sel <- ifelse(hab_depth == "depthsurf",
                       "data/env/current/thetao_baseline_depthsurf_mean.tif",
                       "data/env/current/thetao_baseline_depthmean_mean.tif")
    
    current <- rast(temp_sel)
    
    envelope <- get_envelope(sp_sel_data, current)
    
    scenarios <- data.frame(
      scenario = c("current", rep(c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585"),
                                  each = 2)),
      year = c(NA, rep(c("dec50", "dec100"), 5))
    )

    # Load SST layers
    sst_layers <- rep(NA, length(scenarios))
    for (sc in seq_len(nrow(scenarios))) {
        if (scenarios$scenario[sc] != "current") {
            pred_layer <- gsub("current", paste0("future/", scenarios$scenario[sc]), temp_sel)
            pred_layer <- gsub("baseline", scenarios$scenario[sc], pred_layer)
            pred_layer <- gsub("_mean", paste0("_", scenarios$year[sc], "_mean"), pred_layer)
        } else {
            pred_layer <- temp_sel
        }
        sst_layers[sc] <- pred_layer
    }

    sst_layers <- rast(sst_layers)

    cli::cat_line(cli::col_cyan("Masking SST layers"))
    masked_data <- get_masked(envelope, sst_layers)
    
    limits <- vector(mode = "list", length = nrow(scenarios))
    
    for (i in seq_len(nrow(scenarios))) {
      
      cli::cat_line(cli::col_cyan(paste("Processing scenario", i, "out of", nrow(scenarios))))
      
      limits_t <- terra::extract(sst_layers[[i]], sp_sel_data[,1:2], ID = F)
      limits[[i]] <- data.frame(
        q_0.05 = round(quantile(limits_t[,1], 0.05), 1),
        q_0.5 = round(quantile(limits_t[,1], 0.5), 1),
        q_0.95 = round(quantile(limits_t[,1], 0.95), 1),
        mean_v = round(mean(limits_t[,1]), 1),
        sd_v = round(sd(limits_t[,1]), 1),
        scenario = scenarios[i,1],
        period = scenarios[i,2]
      )
      rownames(limits[[i]]) <- NULL
      
    }
    
    limits <- do.call("rbind", limits)

    masked_data <- as.int(masked_data)
    
    names(masked_data) <- gsub("_NA", "", paste0(scenarios$scenario, "_", scenarios$year))
    
    if (output_format == "multiband") {
      areas <- terra::expanse(masked_data, unit = "km")
      areas <- data.frame(area = areas[,2])

      writeRaster(masked_data, outfile, datatype = "INT1U")
      obissdm::cogeo_optim(outfile)
    } else if (output_format == "parquet") {
      masked_data_pol <- lapply(seq_len(nlyr(masked_data)), function(x) aggregate(as.polygons(masked_data[[x]], values = F)))
      names(masked_data_pol) <- names(masked_data)
      masked_data_pol <- vect(masked_data_pol)
      
      areas <- lapply(seq_len(length(masked_data_pol)), function(x){
        sel_pol <- masked_data_pol[x]
        ex <- terra::expanse(sel_pol, unit = "km")
        data.frame(area = ex)
      })

      masked_data_pol <- terra::simplifyGeom(masked_data_pol)
    
      masked_data_pol_sf <- sf::st_as_sf(masked_data_pol)
      
      sfarrow::st_write_parquet(masked_data_pol_sf, outfile)

      areas <- do.call("rbind", areas)
    }

    areas$scenario <- scenarios$scenario
    areas$year <- scenarios$year
    
    outfile_json <- gsub("thermenvelope.parquet|thermenvelope.tif", "thermmetrics.json", outfile)
    outfile_json <- gsub("predictions", "metrics", outfile_json)
    
    jsonlite::write_json(list(
      areas = areas,
      limits = limits,
      sst_depth = hab_depth
    ), outfile_json, pretty = T)
    
    return("saved")
    
  } else {
    return("low_number")
  }
  
}



# Fit thermal envelopes for all species ----
# Run in parallel
parallel_workers <- 40
plan(multisession, workers = parallel_workers)
mem_frac <- 0.9 / parallel_workers

results <- future_map(tg_species, function(sp) {
  try(get_thermrange(sp, target_folder = results_folder,
                     model_acro = model_acro,
                     output_format = output_format,
                     mem_frac = mem_frac))
}, .progress = T)

print(results[1:4])
