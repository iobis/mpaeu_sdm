############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
################## Run multiple species distribution models ####################


# Load packages ----
library(obissdm)
library(furrr)
library(progressr)
library(storr)
library(polars)
library(arrow)
library(dplyr)
library(terra)
source("functions/model_species.R")
source("functions/auxiliary_modelfit.R")
set.seed(2023)
handlers("cli")
options("progressr.enable" = TRUE)


# Define settings ----
# This file can also be run as "source". In that case, the default objects
# are used. To enable also a more flexible run, for testing purposes, three 
# arguments can be override by creating the objects on the environment before. 
# `Those are outfolder`, `outacro`, and `sel_species`. The latter have as 
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
outacro <- check_exists("outacro", "wwf")
# What will be modelled?
# 'all' will model all species. Supply a vector of AphiaIDs to filter
sel_species <- check_exists("sel_species", "all") 
# The path for the species data dataset
species_dataset <- "data/species/"
# The path for the species list
species_list <- "data/all_splist_20240319.csv"
# Run in parallel? For avoiding parallel, change both to FALSE
run_parallel <- ifelse(length(sel_species) > 1 | sel_species == "all", TRUE, FALSE)[1]
# Number of cores for parallel processing
n_cores <- 4

# Modelling
# Algorithms to be used
algos <- c("maxent", "rf", "brt", "lasso", "xgboost")
# Should areas be masked by the species depth?
limit_by_depth <- TRUE
# A buffer to be applied to the depth limitation
depth_buffer <- 500
# Assess spatial bias?
assess_bias <- FALSE
# Try to correct spatial bias?
correct_bias <- FALSE

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
cli::cli_inform("Running SDM for {length(species_list)} species")
cli::cli_inform("Chosen algorithm{?s}: {.val {algos}}")
if (!limit_by_depth) {
  cli::cli_inform(c("x" = "Not limiting by depth"))
} else {
  cli::cli_inform(c("v" = "Limiting by depth using {.val {depth_buffer}} buffer"))
}
cli::cat_line()
cli::cli_inform("Outputing to folder {.path {outfolder}} using acronym {.val {outacro}}")
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
pmod <- function(sp, gp, sdat, outf, outac, alg, lmd, lmd_buf, assb, corb, p) {
  
  p()
  
  if (st$exists(sp)) {
    if (st$get(as.character(sp)) %in% c("done", "succeeded")) {
      to_do <- FALSE
    } else {
      to_do <- TRUE
    }
  } else {
    to_do <- TRUE
  }
  
  if (to_do) {
    fit_result <- try(model_species(species = sp,
                                    group = gp,
                                    species_dataset = sdat,
                                    outfolder = outf,
                                    outacro = outac,
                                    algorithms = alg,
                                    limit_by_depth = lmd,
                                    depth_buffer = lmd_buf,
                                    assess_bias = assb,
                                    correct_bias = corb),
                      silent = T)
    
    if (!inherits(fit_result, "try-error")) {
      st$set(sp, fit_result)
    } else {
      st$set(sp, "failed")
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
                          sdat = species_dataset, outf = outfolder,
                          outac = outacro, alg = algos, lmd = limit_by_depth,
                          lmd_buf = depth_buffer, assb = assess_bias,
                          corb = correct_bias,
                          p = p, .options = furrr_options(seed = T))
  })
} else {
  with_progress({
    p <- progressor(steps = nrow(species_list))
    result <- purrr::map2(species_list$taxonID, species_list$sdm_group, pmod,
                          sdat = species_dataset, outf = outfolder,
                          outac = outacro, alg = algos, lmd = limit_by_depth,
                          lmd_buf = depth_buffer, assb = assess_bias,
                          corb = correct_bias,
                          p = p)
  })
}



# Check results ----
# Check if everything was processed
cli::cli_alert_warning("{.val {length(st$list())}} out of {.val {nrow(species_list)}} model{?s} processed.")

# And if so, destroy storr object