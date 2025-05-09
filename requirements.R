############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
############################ Project Requirements ##############################

# Needed packages on CRAN
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
  "sfarrow"
)

# Needed packages on GitHub
git_packages <- c("flexsdm", "prg", "obissdm", "obistools", "biooracler")
git_packages_source <- c(
  "sjevelazco/flexsdm@HEAD",
  "meeliskull/prg/R_package/prg",
  "iobis/mpaeu_msdm",
  "iobis/obistools",
  "bio-oracle/biooracler"
)

# Bioconductor packages
bioc_packages <- c("Rarr")

# Packages with special installation
special_packages <- "polars"

# Create a function to check if is installed
is_package_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Check which ones are not installed and install if needed:
for (i in 1:length(req_packages)) {
  if (!is_package_installed(req_packages[i])) {
    install.packages(req_packages[i])
  }
}

# Check github packages
for (i in 1:length(git_packages)) {
  if (!is_package_installed(git_packages[i])) {
    devtools::install_github(git_packages_source[i])
  }
}

# Check bioconductor packages
for (i in 1:length(bioc_packages)) {
  if (!is_package_installed(bioc_packages[i])) {
    BiocManager::install(bioc_packages[i])
  }
}

is_polar_inst <- is_package_installed(special_packages)
if (!is_polar_inst) {
  Sys.setenv(NOT_CRAN = "true")
  install.packages("polars", repos = "https://rpolars.r-universe.dev")
}

# Install minioclient
minioclient::install_mc(force = TRUE)

# Install Python packages
reticulate::py_install(c("rio-cogeo", "xarray", "zarr"), pip = TRUE)