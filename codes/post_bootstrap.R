############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
##################### SDM - bootstrap results of models ########################

# Load packages -----
library(obissdm)
library(terra)
library(storr)
library(progressr)
set.seed(2023)
source("functions/bootstrap.R")
handlers("cli")
options("progressr.enable" = TRUE)
if (interactive()) {
    cat("\014")
} else {
    system("clear")
}
start_time <- Sys.time()

# Model acronym
acro <- "mpaeu"

# Set species to run
# To run specific species, uncomment below and comment the others
sel_species <- c(124287)
# Otherwise, use species that were effectively modeled
# st_sdm <- storr_rds(paste0(acro, "_storr"))
# st_sdm_done <- st_sdm$list()
# st_sdm_status <- st_sdm$mget(st_sdm_done)
# st_sdm_done <- st_sdm_done[unlist(lapply(st_sdm_status, function(x) x[[1]])) == "succeeded"]
# sel_species <- as.numeric(st_sdm_done)

# Run in parallel? For avoiding parallel, change both to FALSE
run_parallel <- ifelse(length(sel_species) > 1 | sel_species == "all", TRUE, FALSE)[1]
# Number of cores for parallel processing
n_cores <- 4
# Bootstrap target ("all", "best" or a particular algorithm)
boot_mode <- "all"

# Create storr to hold results
st <- storr_rds(paste0(acro, "_boot_storr"))
# If does not exist, add start date
if (!st$exists("startdate")) {
    st$set("startdate", format(Sys.Date(), "%Y%m%d"))
}

# Output for control
cli::cat_rule()
cli::cat_line(cli::col_cyan("MPA EUROPE PROJECT - WP3 Species and biogenic habitat distributions"))
cli::cat_line("Bootstrap distribution models of multiple species in parallel")
cli::cat_line()
cli::cli_inform("Bootstrapping SDMs for {length(sel_species)} species")
cli::cat_line()
if (run_parallel) {
    cli::cli_alert_warning("Running in parallel using {n_cores} cores")
} else {
    cli::cli_alert_warning("Running in sequence")
}
cli::cat_line()
cli::cli_alert_info("Using mode {.val {boot_mode}}")
cli::cat_rule()

pmod <- function(species,
                 boot_mode,
                 p) {
    p()

    if (st$exists(species)) {
        if (st$get(as.character(species))[[1]] %in% c("done", "failed")) {
            to_do <- FALSE
        } else {
            to_do <- TRUE
        }
    } else {
        to_do <- TRUE
    }

    if (to_do) {
        st$set(species, "running")

        fit_result <- try(
            bootstrap_sp(species = species, target = boot_mode),
            silent = T
        )

        if (!inherits(fit_result, "try-error")) {
            st$set(species, fit_result)
        } else {
            st$set(species, list(
                status = "failed",
                error = fit_result
            ))
        }
    }

    return(invisible(NULL))
}

# Run models according to the strategy
if (run_parallel) {
    plan(multisession, workers = n_cores)

    with_progress({
        p <- progressor(steps = length(sel_species))
        result <- future_map(sel_species, pmod,
            boot_mode = boot_mode,
            p = p, .options = furrr_options(seed = T)
        )
    })
} else {
    with_progress({
        p <- progressor(steps = length(sel_species))
        result <- purrr::map(sel_species, pmod,
            boot_mode = boot_mode,
            p = p
        )
    })
}

# END
