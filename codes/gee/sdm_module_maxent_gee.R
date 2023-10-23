# Load packages/functions ----
library(rgee)
library(sf)
library(terra)
library(obissdm)
source("functions/load_env.R")
set.seed(2023)


# Set configuration ----
# Authenticate to Google Earth Engine
ee_Authenticate()

# Initialize the Google Earth Engine session
ee_Initialize()


# Prepare data locally ----
# List and load environmental layers
env_vars <- c("thetao-mean")

curr <- load_env(env_vars)

# Load study area shapefile
starea <- vect("data/shapefiles/mpa_europe_starea_v2.shp")
exp_starea <-  ext(-41, 47, 20, 89) # ext(starea) +- 5 degrees

# Crop to the expanded area (only the one that will be used for
# sampling the background)
curr <- crop(curr, exp_starea)

# Remove Red sea
mregions <- mregions::mr_shp("MarineRegions:iho")
mregions <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

curr <- mask(curr, mregions, inverse = T)

# Load species data
key <- 101
pts <- arrow::read_parquet("data/virtual_species/key=101/occurrences_high_rep1.parquet")

# Create background points and generate blocks for cross-validation
block_list <- list(
  spatial_grid = rast(ext(curr), resolution = 5),
  spatial_lat = rast(ext(curr), ncol = 1, nrow = 30)
)

sp_data <- mp_prepare_data(pts, species_id = paste0("species", key),
                           env_layers = curr,
                           quad_number = 150000
)

sp_data <- mp_prepare_blocks(sp_data, method = "manual", manual_shp = block_list,
                             n_iterate = 300, verbose = F)

# Convert data to format for use in Google Earth Engine
sp_data_ee <- cbind(sp_data$coord_training,
                    presence = sp_data$training$presence,
                    spatial_grid = sp_data$blocks$folds$spatial_grid,
                    spatial_lat = sp_data$blocks$folds$spatial_lat,
                    random = sp_data$blocks$folds$random)

sp_data_ee <- st_as_sf(sp_data_ee, coords = c("decimalLongitude", "decimalLatitude"),
                       crs = crs(curr))


# Get best spatial grid blocks ----
n_iterations <- 100
k_folds <- 5

folds <- lapply(1:n_iterations,
                function(x){sample(1:k, nrow(spatial_grid), replace = T)})

folds_ee <- ee$List(folds)

it_ind <- ee$List(0:(n_iterations - 1))

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
  
  # Join and get only information
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


# Sample background points ----
region_rast <- ee$Image("projects/ee-silascprincipe/assets/env_data")

zones <- region_rast$select("b1")$gt(0)
zones <- zones$updateMask(zones$neq(0))

region_pol <- zones$addBands(region_rast$select("b1"))$reduceToVectors(
  crs = region_rast$projection(),
  geometryType = 'polygon',
  reducer = ee$Reducer$mean()
)

# Map$addLayer(region_pol)
back_points <- ee$FeatureCollection$randomPoints(
  region = region_pol$geometry(),
  points = 1000,
  seed = 2,
  maxError = 1)

back_points <- zones$sample(
  region = region_pol$geometry(),
  #scale = 5565.975,
  numPixels = 1000,
  seed = 2,
  geometries = TRUE)

back_points <- zones$sample(
  region = zones$geometry(),
  geometries = TRUE
)

back_points <- back_points$randomColumn(seed = 2)

perc_to_get <- ee$Number((100000 * 100))$divide(ee$Number(back_points$size()))$ceil()

back_points_b <- back_points$filter(ee$Filter$lte("random", ee$Number(perc_to_get$divide(100))))

ee_print(back_points_b)

points_id <- back_points$select("b1")$toList(back_points$size())

new_ind <- ee$List$sequence(0, ee$Number(back_points$size())$subtract(1))

back_points <- back_points$remap(
  points_id, new_ind, "index"
)

ee_print(back_points)

Map$addLayer(back_points)

# Get folds index

# Merge background and points ----

# Extract information at the points ----
sp_data <- env_layers$sampleRegions(collection = sp_points)

# Get best tune ----
# Set parameters
tune_params <- expand.grid(remult = seq(0.5, 2, 0.5), features = c("lq", "lqh"),
                           stringsAsFactors = F)

tune_params <- unique(tune_params)

tune_params$features <- ifelse(grepl("h", tune_params$features), TRUE, FALSE)

tune_params <- lapply(1:nrow(tune_params), function(x){
  regu <- tune_params[x, 1]
  hing <- tune_params[x, 2]
  return(list(regu, hing))
})

# Convert to ee format
tune_list <- ee$List(tune_params)

# Creates a function to get AUC
get_auc_tun <- function(featurecol) {
  
  predicted <- ee$FeatureCollection(featurecol)$aggregate_array("probability")
  response <- ee$FeatureCollection(featurecol)$aggregate_array("presence")
  
  thresholds <- ee$List$sequence(0, 1, 0.01)
  
  thresh_res <- function(th) {
    pred_class <- ee$Array(predicted)$gte(th)
    actu_class <- ee$Array(response)#$gte(th)
    ee$List(list(pred_class, actu_class))
  }
  
  thresholded <- thresholds$map(rgee::ee_utils_pyfunc(thresh_res))
  
  get_classind <- function(classified) {
    
    predvals <- ee$Array(ee$List(classified)$get(0))
    truevals <- ee$Array(ee$List(classified)$get(1))
    
    newvar <- ee$Array(truevals$multiply(2)$add(predvals))

    tp <- newvar$eq(3)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
    fp <- newvar$eq(1)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
    tn <- newvar$eq(0)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))
    fn <- newvar$eq(2)$reduce(reducer = ee$Reducer$sum(), axes = ee$List(list(0L)))

    tpr <- tp$divide(tp$add(fn)) # TP/TP+FN / Recall
    tnr <- tn$divide(tn$add(fp)) # TN/TN+FP
    fpr <- fp$divide(fp$add(tn)) # FP/FP+TN
    pre <- tp$divide(tp$add(fp)) # TP/TP+FP

    ee$Dictionary(list(tp = tp, fp = fp, tn = tn, fn = fn, tpr = tpr, tnr = tnr,
                              fpr = fpr, pre = pre))
    
  }
  
  classif <- thresholded$map(rgee::ee_utils_pyfunc(get_classind))
  
  FPR <- classif$map(rgee::ee_utils_pyfunc(function(dict){
    value = ee$Dictionary(dict)$get('fpr')
    ee$Array(value)$get(ee$List(list(0)))
  }))
  
  TPR <- classif$map(rgee::ee_utils_pyfunc(function(dict){
    value = ee$Dictionary(dict)$get('tpr')
    ee$Array(value)$get(ee$List(list(0)))
  }))
  
  # Concatenate FPR and TPR arrays along axis 1
  X <- ee$Array(FPR) # tentar nao usar array, ver amanha
  Y <- ee$Array(TPR)
  X1 <- X$slice(0,1)$subtract(X$slice(0,0,-1))
  Y1 <- Y$slice(0,1)$add(Y$slice(0,0,-1))
  auc <- X1$multiply(Y1)$multiply(0.5)$reduce('sum', ee$List(list(0)))$abs()$toList()$get(0)
  
  #auc_feat <- ee$Feature(NULL, list(auc = auc))
  
  #return(ee$Number(auc))
  ee$Feature(NULL, list(auc = auc))
  
}

# Creates function to perform tuning
tune_model <- function(listval) {
  
  # Get parameters
  all_parameters <- ee$List(listval)
  
  regmult <- all_parameters$get(0)
  hingefeat <- all_parameters$get(1)
  
  cv_model <- function(fold) {
    
    training = training_dat$filter(ee$String(ee$String("fold != ")$cat(ee$Number(fold)$format())))
    testing = training_dat$filter(ee$String(ee$String("fold == ")$cat(ee$Number(fold)$format())))
    
    maxclassifier <- ee$Classifier$amnhMaxent(
      autoFeature = FALSE,
      linear = TRUE,
      quadratic = TRUE,
      product = FALSE,
      threshold = FALSE,
      hinge = hingefeat, # Only hinge variable
      addSamplesToBackground = FALSE,
      betaMultiplier = ee$Number(regmult) # regularization multiplier
    )$train(
      features = training,
      classProperty = 'presence',
      inputProperties = image$bandNames()
    )
    
    testpred <- ee$FeatureCollection(testing$classify(maxclassifier))
    
    auc_val <- get_auc_tun(testpred)
    
    return(auc_val)
    
  }
  
  fold_list <- ee$List$sequence(1, 2, 1)
  
  auc_cv <- fold_list$map(rgee::ee_utils_pyfunc(cv_model))
  
  mean_auc <- ee$FeatureCollection(auc_cv)$
    aggregate_mean("auc")
  
  return(ee$Feature(NULL, list(auc = mean_auc)))
}

tune_results <- tune_list$map(rgee::ee_utils_pyfunc(tune_model))

tune_results_fc <- ee$FeatureCollection(tune_results)

tr_list <- ee$List(tune_results_fc$aggregate_array("auc"))

assign_list <- ee$List$sequence(0, 7, 1)

tune_res_fin <- assign_list$map(
  rgee::ee_utils_pyfunc(function(index){
    sel_list <- ee$List(tune_list$get(index))
    regmult <- sel_list$get(0)
    hingefeat <- sel_list$get(1)
    auc_val <- tr_list$get(index)
    return(ee$Feature(NULL, list(auc = auc_val,
                                 regmult = regmult,
                                 hingefeat = hingefeat)))
  })
)

best_tune <- ee$FeatureCollection(tune_res_fin)$
  filter(ee$Filter$eq("auc", ee$FeatureCollection(tune_res_fin)$
                        aggregate_min("auc")))


# Cross validate best model ----
regmult <- ee$List(best_tune$aggregate_array("regmult"))$get(0)
hingefeat <- ee$List(best_tune$aggregate_array("hingefeat"))$get(0)

cv_model_final <- function(fold) {
  
  training = training_dat$filter(ee$String(ee$String("fold != ")$cat(ee$Number(fold)$format())))
  testing = training_dat$filter(ee$String(ee$String("fold == ")$cat(ee$Number(fold)$format())))
  
  maxclassifier <- ee$Classifier$amnhMaxent(
    autoFeature = FALSE,
    linear = TRUE,
    quadratic = TRUE,
    product = FALSE,
    threshold = FALSE,
    hinge = hingefeat, # Only hinge variable
    addSamplesToBackground = FALSE,
    betaMultiplier = ee$Number(regmult) # regularization multiplier
  )$train(
    features = training,
    classProperty = 'presence',
    inputProperties = image$bandNames()
  )
  
  testpred <- ee$FeatureCollection(testing$classify(maxclassifier))
  
  # testpred_sel <- testpred$select(
  #   propertySelectors = list("fold", "presence", "probability"),
  #   retainGeometry = FALSE
  # )
  testpred_sel <- testpred$select(
    list("fold", "presence", "probability")
  )
  
  return(ee$FeatureCollection(testpred_sel))
  
}

fold_list <- ee$List$sequence(1, 2, 1)

final_m_cv <- fold_list$map(rgee::ee_utils_pyfunc(cv_model_final))

final_m_cv <- ee$FeatureCollection(final_m_cv)$flatten()


# Train full model ----
maxclassifier_final <- ee$Classifier$amnhMaxent(
  autoFeature = FALSE,
  linear = TRUE,
  quadratic = TRUE,
  product = FALSE,
  threshold = FALSE,
  hinge = hingefeat, # Only hinge variable
  addSamplesToBackground = FALSE,
  betaMultiplier = ee$Number(regmult) # regularization multiplier
)$train(
  features = training_dat,
  classProperty = 'presence',
  inputProperties = image$bandNames()
)

# Predict to full data ----
full_data_pred <- ee$FeatureCollection(training_dat$classify(maxclassifier))

full_data_pred <- full_data_pred$select(
  list("fold", "presence", "probability")
)


# Predict to scenarios ----
for (i in 1:length(scenarios)) {
  
}


# Export predictions ----
for (i in 1:length(scenarios)) {
  
}


# Save objects for further processing ----
# Save best tune

# Save cross-validation metrics

# Save full model metrics 


# Get final CV metrics ----











# DRAFT CODES ----
# 
# best_auc <- ee$FeatureCollection(auc_cv)$
#   filter(ee$Filter$eq("auc", ee$FeatureCollection(auc_cv)$
#                         aggregate_min("auc")))
# 
# auc_array <- ee$FeatureCollection(auc_cv)$aggregate_array("auc")
# 
# for (i in 1:3) {
#   auc_array <- ee$List(list(auc_array, ee$Array(auc_cv$get(1))))
# }
# 
# ee$Array(auc_array)$getInfo()
# 
# first_val <- ee$List(list(auc_cv$get(0)))
# 
# iter_function <- rgee::ee_utils_pyfunc(
#   function (current_element, previous_result){
#     value = ee$List(current_element)
#     ee$List(previous_result)$add(value)
#   }
# )
# 
# resulting_list <- auc_cv$iterate(iter_function, first_val)
# resulting_list_f <- ee$List(ee$List(resulting_list)$slice(1,ee$List(resulting_list)$size()))
# 
# ee$Array(list(resulting_list_f))$getInfo()
# 
# var resulting_array = ee.Array(resulting_list);
# print('resulting_array', resulting_array);
# 
# 
# 
# 
# imageClassified <- image$classify(maxclassifier)
# 
# return(ee$List(imageClassified))
# 
# 
# 
# Map$centerObject(image, 9)
# a<- Map$addLayer(
#   ee$Image(teste$get(0)),
#   visParams = list(bands = 'probability', min = 0, max = 1),
#   name = "Probability b"
# )
# b<- Map$addLayer(
#   ee$Image(teste$get(4)),
#   visParams = list(bands = 'probability', min = 0, max = 1),
#   name = "Probability b"
# )
# a|b
# 
# 
# 
# trainingData = ee$FeatureCollection(list(
#   ee$Feature(ee$Geometry$Point(c(-122.39567, 38.02740)), list(presence = 1, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.68560, 37.83690)), list(presence = 1, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.59755, 37.92402)), list(presence = 0, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.47137, 37.99291)), list(presence = 0, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.52905, 37.85642)), list(presence = 0, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.03010, 37.66660)), list(presence = 0, fold = 1)),
#   ee$Feature(ee$Geometry$Point(c(-122.99567, 38.02740)), list(presence = 1, fold = 2)),
#   ee$Feature(ee$Geometry$Point(c(-122.98560, 37.83690)), list(presence = 1, fold = 2)),
#   ee$Feature(ee$Geometry$Point(c(-122.92905, 37.85642)), list(presence = 0, fold = 2)),
#   ee$Feature(ee$Geometry$Point(c(-122.47137, 37.99291)), list(presence = 0, fold = 2)),
#   ee$Feature(ee$Geometry$Point(c(-122.52905, 37.85642)), list(presence = 0, fold = 2)),
#   ee$Feature(ee$Geometry$Point(c(-122.93010, 37.66660)), list(presence = 0, fold = 2))
# ))
# 
# image = ee$Image('LANDSAT/LC08/C02/T1_L2/LC08_044034_20200606')$
# select(list('.._B.*'))
# 
# training_dat = image$sampleRegions(collection = trainingData, scale= 30)
# 
# classifier = ee$Classifier$amnhMaxent(
#   autoFeature = FALSE,
#   linear = TRUE,
#   quadratic = FALSE,
#   product = FALSE,
#   threshold = FALSE,
#   hinge = FALSE,
#   addSamplesToBackground = FALSE
# )$train(
#   features = training,
#   classProperty = 'presence',
#   inputProperties = image$bandNames()
# )
# 
# imageClassified = image$classify(classifier)
# 
# 
# Map$centerObject(image, 9)
# map_a <- Map$addLayer(
#   imageClassified_a,
#   visParams = list(bands = 'probability', min = 0, max = 1),
#   name = "Probability a"
# )
# map_b <- Map$addLayer(
#   imageClassified_b,
#   visParams = list(bands = 'probability', min = 0, max = 1),
#   name = "Probability b"
# )
# map_a | map_b
# 
# 
# 
# 
# sp_data <- env_layers$sampleRegions(collection = sp_points)
# 
# sp_data <- sp_data$randomColumn(
#   distribution = "uniform"
# )
# 
# # Split out ~80% of the sample for training the classifier.
# 
# cv_model <- function(fold) {
#   
#   training <- sp_data$filter()
#   
#   max_class <- ee$Classifier$amnhMaxent(
#     autoFeature = FALSE,
#     linear = TRUE,
#     quadratic = TRUE,
#     product = FALSE,
#     threshold = FALSE,
#     hinge = FALSE,
#     addSamplesToBackground = FALSE})$
#     train(
#       features = training,
#       classProperty = 'presence',
#       inputProperties = env_layers.bandNames()
#     )
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# training <- sp_data$filter("random < 0.8")
# testing <- sp_data$filter("random > 0.8")
# 
# 
# 
# 
# sampleFeatures <- ee$FeatureCollection(
#   list(ee$Feature(NULL, list(predicted = 0.2, actual = 0)),
#   ee$Feature(NULL, list(predicted = 0.7, actual = 1)),
#   ee$Feature(NULL, list(predicted = 0.8, actual = 1)),
#   ee$Feature(NULL, list(predicted = 0.3, actual = 0)),
#   ee$Feature(NULL, list(predicted = 0.9, actual = 1)))
# )
# 
# 
# get_auc <- function(featurecol) {
#   
#   predicted <- featurecol$aggregate_array("predicted")
#   response <- featurecol$aggregate_array("actual")
#   
#   thresholds <- ee$List$sequence(0, 1, 0.01)
#   
#   thresh_res <- function(th) {
#     pred_class <- ee$Array(predicted)$gte(th)
#     actu_class <- ee$Array(actual)#$gte(th)
#     ee$List(pred_class, actu_class)
#   }
#   
#   thresholded <- thresholds$map(thresh_res)
#   
#   get_classind <- function(classified) {
#     
#     predvals <- ee$Array(ee$List(classified).get(0))
#     truevals <- ee$Array(ee$List(classified).get(1))
#     
#     newvar <- truevals.multiply(2).add(predvals);
#     
#     tp <- newvar$eq(3)$reduce(reducer = ee$Reducer$sum(), axes = 0L)
#     fp <- newvar$eq(1)$reduce(reducer = ee$Reducer$sum(), axes = 0L)
#     tn <- newvar$eq(0)$reduce(reducer = ee$Reducer$sum(), axes = 0L)
#     fn <- newvar$eq(2)$reduce(reducer = ee$Reducer$sum(), axes = 0L)
#     
#     tpr <- tp$divide(tp$add(fn)) # TP/TP+FN / Recall
#     tnr <- tn$divide(tn$add(fp)) # TN/TN+FP
#     fpr <- fp$divide(fp$add(tn)) # FP/FP+TN
#     pre <- tp$divide(tp$add(fp)) # TP/TP+FP
#     
#     return ee.Dictionary(list(tp = tp, fp = fp, tn = tn, fn = fn, tpr = tpr, tnr = tnr,
#       fpr = fpr, pre = pre))
#     
#   }
#   
#   classif <- thresholded$map(get_classind)
#   
#   fpr_array <- ee$Array
#   
#   
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###
# var sampleFeatures = ee.FeatureCollection([
#   ee.Feature(null, {predicted: 0.2, actual: 0}),
#   ee.Feature(null, {predicted: 0.7, actual: 1}),
#   ee.Feature(null, {predicted: 0.8, actual: 1}),
#   ee.Feature(null, {predicted: 0.3, actual: 0}),
#   ee.Feature(null, {predicted: 0.9, actual: 1})
# ]);
# 
# var predicted = sampleFeatures.aggregate_array('predicted');
# var actual = sampleFeatures.aggregate_array('actual');
# 
# var thresholds = ee.List.sequence(0, 1, 0.01);
# 
# print(predicted)
# 
# var rocPoints = function(th){
#   var pred_class = ee.Array(predicted).gte(th);
#   var actu_class = ee.Array(actual).gte(th);
#   return [pred_class, actu_class];
# };
# 
# var thresholded = thresholds.map(rocPoints);
# 
# print(thresholded.get(10))
# 
# 
# /*
#   var getClass = function(classified) {
#     var pred = ee.List(classified).get(0);
#     var actu = ee.List(classified).get(1);
#     var extract = function(actual){
#       return(ee.Array(pred).eq(actual));
#     };
#     var comparison = extract(actu);
#     return comparison;
#   };
# */
#   
#   var getClass = function(classified){
#     var predictedValues = ee.List(classified).get(0);
#     var trueValues = ee.List(classified).get(1);
#     
#     var flatTrue = ee.Array(trueValues);
#     var flatPredicted = ee.Array(predictedValues);
#     
#     var newvar = flatTrue.multiply(2).add(flatPredicted);
#     
#     var tp = newvar.eq(3).reduce({reducer: ee.Reducer.sum(), axes: [0]});
#     var fp = newvar.eq(1).reduce({reducer: ee.Reducer.sum(), axes: [0]});
#     var tn = newvar.eq(0).reduce({reducer: ee.Reducer.sum(), axes: [0]});
#     var fn = newvar.eq(2).reduce({reducer: ee.Reducer.sum(), axes: [0]});
#     
#     //return ee.Dictionary({tp: tp, fp: fp, tn: tn, fn: fn})
#     
#     var tpr = tp.divide(tp.add(fn)); //TP/TP+FN
#     var tnr = tn.divide(tn.add(fp)); //TN/TN+FP
#     
#     return ee.Dictionary({tp: tp, fp: fp, tn: tn, fn: fn, tpr: tpr, tnr: tnr});
#   };
# 
# var classif = thresholded.map(getClass);
# 
# print(classif);
# 
# print(ee.Dictionary(classif.get(10)).get('tpr'))
# 
# 
# var FPR = classif.map(function(dict){
#   var value = ee.Dictionary(dict).get('tpr');
#   return ee.Array(value).get([0]);
# })
# var TPR = classif.map(function(dict){
#   var value = ee.Dictionary(dict).get('tpr');
#   return ee.Array(value).get([0]);
# })
# 
# print(ee.Array(FPR.flatten()))
# 
# 
# // Concatenate FPR and TPR arrays along axis 1
# var X = ee.Array(FPR);
# var Y = ee.Array(TPR); 
# var X1 = X.slice(0,1).subtract(X.slice(0,0,-1));
# var Y1 = Y.slice(0,1).add(Y.slice(0,0,-1));
# var auc = X1.multiply(Y1).multiply(0.5).reduce('sum',[0]).abs();
# 
# print(auc)