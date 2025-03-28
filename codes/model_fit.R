############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################## Run multiple species distribution models ####################


# Load packages ----
cli::cat_line(cli::col_cyan("This is `obissdm` version ", packageVersion("obissdm"),
                            ". Loading packages..."))
library(obissdm)
library(furrr)
library(progressr)
library(storr)
library(polars)
library(arrow)
library(dplyr)
library(terra)
source("functions/model_species.R")
source("functions/components_model_species.R")
source("functions/auxiliary_modelfit.R")
set.seed(2023)
handlers("cli")
options("progressr.enable" = TRUE)
if (interactive()) { cat("\014") } else { system("clear") }
start_time <- Sys.time()


# Define settings ----
# This file can also be run as "source". In that case, the default objects
# are used. To enable also a more flexible run, for testing purposes, three 
# arguments can be override by creating the objects on the environment before. 
# Those are `outfolder`, `outacro`, and `sel_species`. The latter have as 
# default the value "all", what means that all species in the species list will 
# be modeled. Supplying a vector of AphiaID will filter the species based on the
# supplied values.
#
# We create a small function to perform those checks
check_exists <- function(object_name, if_not) {
  if (exists(object_name)) {
    return(eval(parse(text = object_name)))
  } else {
    return(if_not)
  }
}

# General
# The output folder
outfolder <- check_exists("outfolder", "results")
# An unique code for identifying this model run
outacro <- check_exists("outacro", "mpaeu")
# What will be modelled?
# 'all' will model all species. Supply a vector of AphiaIDs to filter
sel_species <- check_exists("sel_species", "all")
# The path for the species data dataset
species_dataset <- "data/species/"
# The path for the species list
species_list <- recent_file("data", "all_splist")
# Run in parallel? For avoiding parallel, change both to FALSE
run_parallel <- ifelse(length(sel_species) > 1 | sel_species == "all", TRUE, FALSE)[1]
# Number of cores for parallel processing
n_cores <- 40
# Maximum memory used by `terra`
max_mem <- (0.9/n_cores)

# Modelling
# Algorithms to be used
algos <- c("maxent", "rf", "xgboost")
# Personalized options
algo_opts <- obissdm::sdm_options()[algos]
algo_opts$maxent$features <- c("lq", "h")
algo_opts$maxent$remult <- seq_len(4)
algo_opts$xgboost$gamma <- c(0, 4)
algo_opts$xgboost$shrinkage <- c(0.1, 0.3)
algo_opts$xgboost$scale_pos_weight <- c("balanced", "equal")
# Should areas be masked by the species depth?
limit_by_depth <- TRUE
# A buffer to be applied to the depth limitation
depth_buffer <- 50
# Assess spatial bias?
assess_bias <- TRUE
# Quadrature size
quad_samp <- 0.01 # 1% of the total number of points
# Target metric
tg_metric <- "cbi"
# Metric threshold
tg_threshold <- 0.3

# Create storr to hold results
st <- storr_rds(paste0(outacro, "_storr"))
# If does not exist, add start date
if (!st$exists("startdate")) {
  st$set("startdate", format(Sys.Date(), "%Y%m%d"))
}
# Should the storr object be destructed at the end if the model succeeded?
destroy_storr <- FALSE


# Define species ----
species_list <- read.csv(species_list)
species_list <- get_listbygroup(species_list)

if (length(sel_species) > 1 | !any("all" %in% sel_species)) {
  species_list <- species_list %>%
    filter(taxonID %in% sel_species)
}


# Output for control
cli::cat_rule()
cli::cat_line(cli::col_cyan("MPA EUROPE PROJECT - WP3 Species and biogenic habitat distributions"))
cli::cat_line("Run species distribution models for multiple species in parallel")
cli::cat_line()
cli::cli_inform("Running SDM for {nrow(species_list)} species")
cli::cli_inform("Chosen algorithm{?s}: {.val {algos}}")
if (!limit_by_depth) {
  cli::cli_inform(c("x" = "Not limiting by depth"))
} else {
  cli::cli_inform(c("v" = "Limiting by depth using {.val {depth_buffer}} buffer"))
}
cli::cat_line()
cli::cli_inform("Outputting to folder {.path {outfolder}} using acronym {.val {outacro}}")
if (dir.exists(outfolder)) {
  cli::cli_alert_warning("Folder already exists, files may be overwritten!")
} else {
  fs::dir_create(outfolder)
}
if (run_parallel) {
  cli::cli_alert_warning("Running in parallel using {n_cores} cores")
} else {
  cli::cli_alert_warning("Running in sequence")
}
cli::cat_rule()




# Fit models ----
# Create a function to save results in a storr object for better control
pmod <- function(species,
                 group,
                 species_dataset,
                 outfolder,
                 outacro,
                 algorithms,
                 algo_opts = NULL,
                 limit_by_depth,
                 depth_buffer,
                 assess_bias,
                 quad_samp,
                 max_mem,
                 p) {
  
  p()
  
  if (st$exists(species)) {
    if (st$get(as.character(species))[[1]] %in% c("done", "succeeded", "low_data",
                                             "failed", "no_good_model")) {
      to_do <- FALSE
    } else {
      to_do <- TRUE
    }
  } else {
    to_do <- TRUE
  }
  
  if (to_do) {
    st$set(species, "running")
    
    fit_result <- try(model_species(
      species = species,
      group = group,
      species_dataset = species_dataset,
      outfolder = outfolder,
      outacro = outacro,
      algorithms = algos,
      algo_opts = algo_opts,
      limit_by_depth = limit_by_depth,
      depth_buffer = depth_buffer,
      assess_bias = assess_bias,
      post_eval = c("sst", "niche", "hyper"),
      tg_metric = tg_metric,
      tg_threshold = tg_threshold,
      quad_samp = quad_samp,
      max_mem = max_mem,
      verbose = FALSE
    ), silent = T)
    
    if (!inherits(fit_result, "try-error")) {
      st$set(species, fit_result)
    } else {
      st$set(species, list(status = "failed",
                      error = fit_result))
    }
  }
  
  return(invisible(NULL))
}

# Run models according to the strategy
if (run_parallel) {
  plan(multisession, workers = n_cores)
  
  with_progress({
    p <- progressor(steps = nrow(species_list))
    result <- future_map2(species_list$taxonID, species_list$sdm_group, pmod,
                          species_dataset = species_dataset,
                          outfolder = outfolder,
                          outacro = outacro,
                          algorithms = algos,
                          algo_opts = algo_opts,
                          limit_by_depth = limit_by_depth,
                          depth_buffer = depth_buffer,
                          assess_bias = assess_bias,
                          quad_samp = quad_samp,
                          max_mem = max_mem,
                          p = p, .options = furrr_options(seed = T))
  })
} else {
  with_progress({
    p <- progressor(steps = nrow(species_list))
    result <- purrr::map2(species_list$taxonID, species_list$sdm_group, pmod,
                          species_dataset = species_dataset,
                          outfolder = outfolder,
                          outacro = outacro,
                          algorithms = algos,
                          algo_opts = algo_opts,
                          limit_by_depth = limit_by_depth,
                          depth_buffer = depth_buffer,
                          assess_bias = assess_bias,
                          quad_samp = quad_samp,
                          max_mem = max_mem,
                          p = p)
  })
}



# Check results ----
# Check if everything was processed
sp_done <- st$list()
sp_done <- length(sp_done[sp_done %in% species_list$taxonID])
cli::cli_alert_warning("{.val {sp_done}} out of {.val {nrow(species_list)}} model{?s} processed.")

# Save session info
fs::dir_create("data/log")
writeLines(c(capture.output(devtools::session_info()),
             capture.output(
              glue::glue(""),
              glue::glue("Model fitting information"),
              glue::glue("Acronym: {outacro}"),
              glue::glue("Start time: {start_time}"),
              glue::glue("End time: {Sys.time()}"),
              glue::glue("Number of species: {nrow(species_list)}"),
              glue::glue("Number of species processed: {sp_done}"),
              glue::glue("obissdm version: {packageVersion('obissdm')}")
             )),
           paste0("data/log/", outacro, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_sessioninfo.txt"))

# And if so, destroy storr object
if (destroy_storr) {
  cli::cli_alert_warning("Destroying `storr` object at {.path {paste0(outacro, '_storr')}}")
  st$destroy()
}