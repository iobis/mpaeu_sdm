#' Load environmental layers downloaded from Bio-ORACLE v3
#'
#' @param ... names of the layers (can also be a vector of names) including the 
#'  variant (min, mean, max, etc) separated by a "-",
#'  e.g. "thethao-mean" will load mean sea surface temperature.
#'  If empty, then terrain_vars should be supplied.
#' @param depth the depth of the layer. If one single value is supplied, it is
#'  used for all variables. If a different depth is intended for each varibale,
#'  then supply a vector of the same length as variables.
#' @param scenario the scenario for which load the variables.
#'  Should be "current" or one of the SSP (e.g. "ssp126")
#' @param decade in case of future scenarios, the decade to be loaded 
#'  (in the format "dec#", being # the decaded). Ignored when loading current 
#'  period or terrain variables.
#' @param terrain_vars to load one or more terrain variable, supply its full name 
#'  (e.g. bathymetry_mean)
#'
#' @return loaded layers (SpatRaster)
#' @export
#'
#' @examples
#' \dontrun{
#' load_env(c("po4-min", "thetao-max"))
#' }
load_env <- function(..., depth = "surf", scenario = "current",
                     decade = "dec100", terrain_vars = NULL) {
  
  variables <- c(...)
  
  to_load <- c()
  
  if (is.null(variables) & is.null(terrain_vars)) {
    stop("At least one variable or a terrain variable should be supplied.")
  }
  
  if (!is.null(variables)) {
    
    if (length(depth) > 1) {
      if (length(depth) != length(variables)) {
        stop("Depth should be a single value or a vector the same size as variables supplied.")
      }
    }
    
    splited <- strsplit(variables, "-")
    
    variables <- unlist(lapply(splited, function(x) x[1]))
    variant <- unlist(lapply(splited, function(x) x[2]))
    
    if (scenario == "current") {
      to_load <- c(to_load,
                   paste0("data/env/current/", variables,
                          "_baseline_depth", depth, "_", variant, ".tif"))
    }
    if (scenario != "current") {
      to_load <- c(to_load,
                   paste0("data/env/future/", scenario, "/", variables,
                          "_", scenario, "_depth", depth,
                          "_", decade, "_",
                          variant, ".tif"))
    }
  }
  
  
  if (!is.null(terrain_vars)) {
    to_load <- c(to_load, paste0("data/env/terrain/", terrain_vars, ".tif"))
  }
  
  return(terra::rast(to_load))
  
}
