# Internal functions - to document

# Convert eco info to habitat depth
hab_to_depth <- function(hab, default = "depthsurf") {
  
  if (grepl("benthic|demersal|bottom", hab)) {
    dstrata <- "depthmax"
  } else if (grepl("pelagic_surface", hab)) {
    dstrata <- "depthsurf"
  } else if (grepl("NOT_FOUND", hab)) {
    dstrata <- default
  } else {
    dstrata <- "depthmean"
  }
  
  return(dstrata)
}

# Split dataset into eval and fit
split_ds <- function(sp_occ,
                          what = "fit",
                          only_coords = TRUE) {
  
  if (grepl("fit", what)) {
    to_get <- "fit_points"
  } else if (grepl("eval", what)) {
    to_get <- "eval_points"
  } else {
    stop("What should be one of `fit_points` or `eval_points`")
  }
  
  pts <- sp_occ %>%
    filter(data_type == to_get) %>%
    as.data.frame()
  
  if (nrow(pts) == 0) {
    pts <- NULL
  } else if (only_coords) {
    pts <- pts[,c("decimalLongitude", "decimalLatitude")]
  }
  
  return(pts)
}

