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

# Settings ----
model_acro <- "inteval"

# Prepare adapted function to get thermal envelope
# This function is based on the `get_thermal_envelope` of the package `speedy`
# but use the data we used for SDMs instead
# https://github.com/iobis/rspeedy

get_envelope <- function (sp_pts, env) {
 
  # Use same standards as speedy
  sf::sf_use_s2(FALSE)
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
  
  if (length(dist) > max_points) {
    points <- terra::vect(sf::st_sample(sf::st_as_sf(dist), max_points))
  } else {
    points <- dist
  }
  
  temperatures <- terra::extract(env, points)
  temperatures <- temperatures[,2]
  
  if (length(temperatures) < min_temperatures) {
    return(NULL)
  }
  
  kd <- ks::kde(temperatures, h = kde_bandwidth)
  percentiles <- ks::qkde(c(0.01, 0.99), kd)
  
  return(percentiles)
}

get_masked <- function(percentiles, layer) {
  
  masked <- layer
  masked[layer >= percentiles[1] & layer <= percentiles[2]] <- 1
  masked[layer < percentiles[1] | layer > percentiles[2]] <- NA
  
  return(masked)
}

# Create a function to generate the thermal range maps
get_thermrange <- function(species, skip_done = TRUE) {
  
  # Load species data
  sp_data <- open_dataset("data/species/")
  sp_code <- species
  
  sp_sel_data <- sp_data %>%
    filter(taxonID == sp_code) %>%
    collect()
  
  outfile <- paste0("results/taxonid=", species, "/model=", model_acro,
                    "/predictions/taxonid=", species, "_model=",
                    model_acro, "_what=thermenvelope.parquet")
  
  if (file.exists(outfile) & skip_done) {
    return("already_saved")
  }
  
  if (nrow(sp_sel_data) > 15) {
    
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
    
    masked_data <- list()
    limits <- list()
    
    for (i in 1:nrow(scenarios)) {
      
      cli::cat_line(cli::col_cyan(paste("Processing scenario", i, "out of", nrow(scenarios))))
      
      if (scenarios$scenario[i] != "current") {
        pred_layer <- gsub("current", paste0("future/", scenarios$scenario[i]), temp_sel)
        pred_layer <- gsub("baseline", scenarios$scenario[i], pred_layer)
        pred_layer <- gsub("_mean", paste0("_", scenarios$year[i], "_mean"), pred_layer)
      } else {
        pred_layer <- temp_sel
      }
      pred_layer <- rast(pred_layer)
      
      masked_data[[i]] <- get_masked(envelope, pred_layer)
      limits_t <- terra::extract(pred_layer, sp_sel_data[,1:2], ID = F)
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
    
    masked_data <- rast(masked_data)
    masked_data <- as.int(masked_data)
    
    names(masked_data) <- gsub("_NA", "", paste0(scenarios$scenario, "_", scenarios$year))
    masked_data_pol <- lapply(1:nlyr(masked_data), function(x) aggregate(as.polygons(masked_data[[x]], values = F)))
    names(masked_data_pol) <- names(masked_data)
    masked_data_pol <- vect(masked_data_pol)
    
    areas <- lapply(1:length(masked_data_pol), function(x){
      sel_pol <- masked_data_pol[x]
      ex <- terra::expanse(sel_pol, unit = "km")
      data.frame(area = ex)
    })
    areas <- do.call("rbind", areas)
    areas$scenario <- scenarios$scenario
    areas$year <- scenarios$year
    
    masked_data_pol <- terra::simplifyGeom(masked_data_pol)
    
    masked_data_pol_sf <- sf::st_as_sf(masked_data_pol)
    
    sfarrow::st_write_parquet(masked_data_pol_sf, outfile)
    
    outfile_json <- gsub("thermenvelope.parquet", "thermmetrics.json", outfile)
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

# Temporary:
tg_species <- c(
  "100803", "100961", "101174", "105883", "107253", "107380", "107392", "107409", "107442", "107455",
  "107518", "111723", "125349", "125366", "126459", "126830", "126868", "127049", "127214", "127380",
  "127393", "130490", "132371", "135146", "135178", "135180", "136290", "137080", "139454", "1398844",
  "140770", "140780", "141436", "141774", "1420120", "145297", "145306", "145307", "145309", "145310",
  "145511", "145540", "145541", "145544", "145548", "145551", "145716", "145723", "145725", "145728",
  "145734", "145794", "145795", "1465385", "146733", "207507", "207516", "210726", "234483", "238025",
  "283704", "286731", "288889", "289939", "371985", "372019", "372502", "383014", "494791", "495082",
  "495309"
)

# Run in parallel
plan(multisession, workers = 4)

results <- future_map(as.numeric(tg_species), function(sp) try(get_thermrange(sp)), .progress = T)

print(results[1:4])