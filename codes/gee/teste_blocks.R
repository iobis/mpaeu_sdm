# Read species data
sp_points <- arrow::read_parquet("data/virtual_species/key=101/occurrences_high_bias_rep1.parquet")
sp_points <- st_as_sf(sp_points[,1:3], coords = c("decimalLongitude", "decimalLatitude"), crs = crs(curr))

# Create an EE object
sp_points_ee <- sf_as_ee(sp_points)

index_list <- seq(0, (nrow(sp_points)-1))

# Random number lists
my_list <- list()

for (i in 1:10) {
  my_list[[i]] <- sample(1:5, nrow(sp_points), replace = T)
}

my_list <- sample(1:5, nrow(sp_points), replace = T)

my_list_ee <- ee$List(my_list)

sp_points_list <- sp_points_ee$toList(sp_points_ee$size())

index_list <- ee$List(index_list)

set_new <- ee_utils_pyfunc(function(id){
  feat <- ee$Feature(sp_points_list$get(id))
  feat$set('idnew', my_list_ee$get(id))
})

new_obj <- index_list$map(set_new)

new_obj_fc <- ee$FeatureCollection(new_obj)

ee_as_sf(new_obj_fc$limit(6))

my_list[1:6]

teste_dic <- list(index = c("a1", "a2"), value = c(0, 2))

teste_dic_ee <- ee$Dictionary(teste_dic)

teste_dic_ee$getInfo()



spatial_grid <- block_list$spatial_grid
spatial_grid[] <- 1:ncell(spatial_grid)
names(spatial_grid) <- "cell"
spatial_grid <- st_as_sf(as.polygons(spatial_grid))
spatial_grid_ee <- sf_as_ee(spatial_grid)

cells_id <- ee$List(spatial_grid$cell)

folds_id <- ee$List(sample(1:5, nrow(spatial_grid), replace = T))

spatial_grid_fold <- spatial_grid_ee$remap(
  cells_id, folds_id, "cell"
)

spat_filter <- ee$Filter$intersects(
  leftField = ".geo",
  rightField = ".geo",
  maxError = 10
)

spat_join <- ee$Join$saveAll(
  matchesKey = "matches"
)

points_info <- spat_join$apply(sp_points, spatial_grid_fold, spat_filter)

match_fun <- ee_utils_pyfunc(function(shape){
  #ee$Feature(shape)$toDictionary()
  ee$Feature(shape)$get("cell")
})

extract_info <- function(point) {
  point <- ee$Feature(point)
  matches <- ee$List(point$get("matches"))
  properties <- matches$map(match_fun)
  point$set("grid_cell", properties)$select(propertySelectors = c("presence", "grid_cell"))
}

points_with_grid <- points_info$map(extract_info)

folds_dist <- points_with_grid$reduceColumns(
  selectors = list("grid_cell", "presence"),
  reducer = ee$Reducer$sum()$group(
    groupField = 0,
    groupName = "fold"
  )
)

folds_metrics <- ee$List(folds_dist$get('groups'))


folds_metrics_final <- folds_metrics$map(ee_utils_pyfunc(
  function(feat){
    ee$Dictionary(feat)$get(key = "sum")
  }
))

new_result <- ee$FeatureCollection(
  ee$Feature(NULL, list(
    combination = ee$Number(1),
    eveness = ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
    all_non_zero = ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero())
  ))
)

ee_as_sf(new_result)

folds_metrics_final$getInfo()

ee_as_sf(points_with_grid$limit(5))




# folds <- list()
# 
# for (i in 1:3) {
#   
#   # Get folds
#   folds[[i]] <- sample(1:5, nrow(spatial_grid), replace = T)
#   
#   # Create folds IDs
#   folds_id <- ee$List(folds[[i]])
#   
#   # Remap the grid
#   spatial_grid_fold <- spatial_grid_ee$remap(
#     cells_id, folds_id, "cell"
#   )
#   
#   # Extract from shapefile
#   # Create a filter
#   spat_filter <- ee$Filter$intersects(
#     leftField = ".geo",
#     rightField = ".geo",
#     maxError = 10
#   )
#   
#   # Create a join
#   spat_join <- ee$Join$saveAll(
#     matchesKey = "matches"
#   )
#   
#   # Join and get only informatio
#   points_info <- spat_join$apply(sp_points, spatial_grid_fold, spat_filter)
# 
#   extract_info <- function(point) {
#     point <- ee$Feature(point)
#     matches <- ee$List(point$get("matches"))
#     properties <- matches$map(ee_utils_pyfunc(function(shape){
#       ee$Feature(shape)$get("cell")
#     }))
#     point$set("grid_cell", properties)$select(propertySelectors = c("presence", "grid_cell"))
#   }
#   points_with_grid <- points_info$map(extract_info)
#   
#   # Calculate distribution
#   folds_dist <- points_with_grid$reduceColumns(
#     selectors = list("grid_cell", "presence"),
#     reducer = ee$Reducer$sum()$group(
#       groupField = 0,
#       groupName = "fold"
#     )
#   )
#   
#   # Extract only relevant field
#   folds_metrics <- ee$List(folds_dist$get('groups'))
#   
#   folds_metrics_final <- folds_metrics$map(ee_utils_pyfunc(
#     function(feat){
#       ee$Dictionary(feat)$get(key = "sum")
#     }
#   ))
#   
#   # Save as a new feature collectio or bind
#   if (i == 1) {
#     new_result <- ee$FeatureCollection(
#       ee$Feature(NULL, list(
#         combination = ee$Number(i),
#         eveness = ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
#         all_non_zero = ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero())
#       ))
#     )
#   } else {
#     new_result <- new_result$merge(
#       ee$FeatureCollection(
#         ee$Feature(NULL, list(
#           combination = ee$Number(i),
#           eveness = ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
#           all_non_zero = ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero())
#         ))
#       )
#     )
#   }
# }
# 
# ee_as_sf(new_result)

folds <- lapply(1:100, function(x){sample(1:5, nrow(spatial_grid), replace = T)})

folds_ee <- ee$List(folds)

it_ind <- ee$List(0:2)

find_best_block <- function(index) {
  
  # Create folds IDs
  folds_id <- ee$List(folds_ee$get(index))
  
  # Remap the grid
  spatial_grid_fold <- spatial_grid_ee$remap(
    cells_id, folds_id, "cell"
  )
  
  # Extract from shapefile
  # Create a filter
  spat_filter <- ee$Filter$intersects(
    leftField = ".geo",
    rightField = ".geo",
    maxError = 10
  )
  
  # Create a join
  spat_join <- ee$Join$saveAll(
    matchesKey = "matches"
  )
  
  # Join and get only informatio
  points_info <- spat_join$apply(sp_points, spatial_grid_fold, spat_filter)
  
  extract_info <- function(point) {
    point <- ee$Feature(point)
    matches <- ee$List(point$get("matches"))
    properties <- matches$map(ee_utils_pyfunc(function(shape){
      ee$Feature(shape)$get("cell")
    }))
    point$set("grid_cell", properties)$select(propertySelectors = c("presence", "grid_cell"))
  }
  points_with_grid <- points_info$map(extract_info)
  
  # Calculate distribution
  folds_dist <- points_with_grid$reduceColumns(
    selectors = list("grid_cell", "presence"),
    reducer = ee$Reducer$sum()$group(
      groupField = 0,
      groupName = "fold"
    )
  )
  
  # Extract only relevant field
  folds_metrics <- ee$List(folds_dist$get('groups'))
  
  folds_metrics_final <- folds_metrics$map(ee_utils_pyfunc(
    function(feat){
      ee$Dictionary(feat)$get(key = "sum")
    }
  ))
  
  # Get final metrics
  # new_result <- ee$List(
  #   list(
  #     ee$Number(index),
  #     ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
  #     ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero())
  #   )
  # )
  
  # new_result <- ee$Feature(NULL, list(
  #   combination = ee$Number(index),
  #   eveness = ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
  #   all_non_zero = ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero())
  # ))
  
  new_result <- ee$List(
    list(ee$Number(index),
         ee$List(folds_metrics_final)$reduce(ee$Reducer$stdDev()),
         ee$List(folds_metrics_final)$reduce(ee$Reducer$allNonZero()))
  )
  
  new_result
}

all_blocks <- it_ind$map(ee_utils_pyfunc(find_best_block))

all_blocks_col <- ee$FeatureCollection(all_blocks$map(ee_utils_pyfunc(
  function(val){
    ee$Feature(NULL, ee$Dictionary$fromLists(
      ee$List(list("combination", "eveness", "all_non_zero")),
      val
    ))
  }
)))

all_blocks_col <- all_blocks_col$
  filter("all_non_zero == 1")

best_block <- 
  ee$Number(all_blocks_col$
              filter(ee$Filter$lte("eveness", all_blocks_col$aggregate_min("eveness")))$
              aggregate_array("combination")$
              get(0))

# Get final remap
spatial_grid_final <- spatial_grid_ee$remap(
  cells_id, ee$List(folds_ee$get(best_block)), "cell"
)

ee_as_sf(spatial_grid_final)
