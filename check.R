############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
######################## Check project Requirements ############################

mes <- function(print_text, danger = F, warning = F, incyan = F) {
  if (incyan) {
    cli::cli_alert_info(cli::col_cyan(print_text))
  } else {
    if (danger) {
      cli::cli_alert_danger(print_text)
    } else if (warning) {
      cli::cli_alert_warning(print_text)
    } else {
      cli::cli_alert_success(print_text)
    }
  }
}

dangers <- 0
warnings <- 0

# Check that folder structure is consistent
de <- lapply(c("data", "functions", "codes"), function(x) dir.exists(x))
if (!all(unlist(de))) {
  mes("Main folders not available. Check if running from root directory and if GitHub repo was correctly cloned.", danger = T)
  dangers <- dangers + 1
} else {
  mes("Main folders available")
}

# Check that configuration file is available:
if (!file.exists("sdm_conf.yml")) {
  mes("Configuration file not available.", danger = T)
  dangers <- dangers + 1
} else {
  mes("Configuration file available.")
}

# Check if species lists are available
if (!any(grepl("all_splist", list.files("data/")))) {
  mes("Species list from MPA Europe not available, you will need to get it or run script for species list.", warning = T)
  warnings <- warnings + 1
} else {
  mes("Species list available.")  
}

# Check with required packages were correctly installed
req_packages <- c(
  "furrr",
  "progressr",
  "storr",
  "arrow",
  "dplyr",
  "terra",
  "stars",
  "cli",
  "fs",
  "purrr",
  "devtools",
  "hypervolume",
  "ggplot2",
  "patchwork",
  "tools",
  "utils",
  "sf",
  "ecospat",
  "ks",
  "ade4",
  "spatstat",
  "ragg",
  "stringr",
  "raster",
  "predicts",
  "worrms",
  "rstudioapi",
  "rgbif",
  "yaml",
  "glue",
  "virtualspecies",
  "lubridate",
  "geohashTools",
  "rfishbase",
  "robis",
  "readr",
  "DBI",
  "duckdb",
  "blockCV",
  "rnaturalearth",
  "isotree",
  "reticulate",
  "sp",
  "ecodist",
  "spThin",
  "modEvA",
  "precrec",
  "rlang",
  "maxnet",
  "glmnet",
  "dismo",
  "randomForest",
  "mgcv",
  "xgboost",
  "lightgbm",
  "leaflet",
  "sys",
  "jsonlite",
  "gbifdb",
  "minioclient",
  "BiocManager",
  "flexsdm",
  "prg",
  "obissdm",
  "obistools",
  "biooracler",
  "Rarr",
  "polars"
)

# Create a function to check if is installed
is_package_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Check which ones are not installed and install if needed:
not_installed <- c()
for (i in 1:length(req_packages)) {
  if (!is_package_installed(req_packages[i])) {
    not_installed <- c(not_installed, req_packages[i])
  }
}

if (length(not_installed) > 0) {
  mes("Package{?s} {.pkg {not_installed}} not installed.", danger = T)
  dangers <- dangers + 1
} else {
  mes("All packages installed.")
}

# Check Python dependencies
req_packages <- c("rio_cogeo", "xarray")
not_installed <- c()
for (i in 1:length(req_packages)) {
  if (!reticulate::py_module_available(req_packages[i])) {
    not_installed <- c(not_installed, req_packages[i])
  }
}

if (length(not_installed) > 0) {
  mes("Python module{?s} {.pkg {not_installed}} not installed.", danger = T)
  dangers <- dangers + 1
} else {
  mes("All Python modules installed.")
}

# Check if rio-cogeo is available from command line
rc <- try(system("rio cogeo --help", intern = T, ignore.stdout = T, ignore.stderr = T), silent = T)
if (inherits(rc, "try-error")) {
  mes("rio-cogeo not available, or not on PATH. Follow instructions from {.url https://cogeotiff.github.io/rio-cogeo/} to install it.", danger = T)
  dangers <- dangers + 1
} else {
  mes("rio-cogeo installed and working.")
}

if (dangers > 0 | warnings > 0) {
  mes("Not all requirements completed. There were {dangers} danger message{?s} and {warnings} warning{?s}.", incyan = T)
} else {
  mes("All requirements completed. Ready to proceed.", incyan = T)
}
