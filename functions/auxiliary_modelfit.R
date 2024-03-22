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
# 
# # Generate log object
# gen_log <- function(algos){
#   
#   log_obj <- list(
#     taxonID = NULL,
#     scientificName = NULL,
#     model_date = NULL,
#     model_acro = NULL,
#     model_fit_points = NULL,
#     model_eval_points = NULL,
#     algorithms = algos,
#     algorithms_parameters = eval(parse(text = paste0(
#       "list(", paste0(algos, "= unclass(sdm_options('", algos, "'))", collapse = ","), ")"
#     ))),
#     model_result = eval(parse(text = paste0(
#       "list(", paste0(algos, "=", "NULL", collapse = ","), ")"
#     ))),
#     model_details = list(
#       ecoregions = NULL,
#       ecoregions_included = NULL,
#       limited_by_depth = NULL,
#       depth_buffer = NULL,
#       block_size = NULL
#     ),
#     obissdm_version = as.character(packageVersion("obissdm"))
#   )
#   
#   return(log_obj)
# }
# 
# 
# 
# # 
# # 
# # 
# # TODO!
# extract_data(sdm_data, what = "coords", from = "train", where = NULL, condition = NULL) {
#   
#   if (what == "coords") {
#     if (is.null(where)) {
#       if (from == "train") {
#         sp_data$coord_training
#       } else {
#         sp_data$coord_eval
#       }
#     } else {
#       if (from == "train") {
#         sp_data$coord_training[sp_data$training[,"presence"] == 1]
#       } else {
#         sp_data$coord_eval
#       }
#     }
#   }
#   
# }