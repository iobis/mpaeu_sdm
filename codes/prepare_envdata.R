### Prepare env data
library(terra)
library(obissdm)
set.seed(2023)

# All croping and masking will be done for each layer
# The depth of cut is registered in the model log file

# Here we assess the VIF of the chosen layers
# We do this for each depth/group combination
# The final decision on which remain is, however, done by hand


# Prepare layers
clip_area <- ext(-41, 47, 20, 89)

# mregions <- mregions::mr_shp("MarineRegions:iho")
# mask_area <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]
mask_area <- vect("data/shapefiles/mpa_europe_starea_v2.shp")

groups <- names(get_conf(what = "groups")$groups)
depth_env <- c("depthsurf", "depthmean")

test_grid <- expand.grid(groups, depth_env, stringsAsFactors = F)


# Create a list to hold the results of vif step
vif_step_list <- list()

for (i in 1:nrow(test_grid)) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 2],
                               load_all = T,
                               verbose = TRUE)
  
  nams <- unique(unlist(env_layers$hypothesis))
  nams <- nams[!grepl("wavefetch", nams)]
  
  env_layers <- terra::subset(env_layers$layers, nams)
  
  # env_layers <- crop(env_layers, mask_area)
  # 
  # env_layers <- mask(env_layers, mask_area, inverse = T)
  
  vif_step_list[[i]] <- usdm::vifstep(env_layers, th = 5)
}

# See which ones have collinearity problems
which_excluded <- lapply(1:length(vif_step_list), function(x) {
  exc <- vif_step_list[[x]]@excluded
  if (length(exc) > 0) {
    cat(x, ">>>", test_grid[x, 1], "|", test_grid[x, 2], "\n")
    cat(exc, "\n\n")
    return(x)
  } else {
    return(invisible(NULL))
  }
})

to_check <- unlist(which_excluded)

to_check


# For those we will get the correlation matrix
vif_list <- list()

for (i in to_check) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 2],
                               load_all = T,
                               verbose = FALSE)
  
  nams <- unique(unlist(env_layers$hypothesis))
  nams <- nams[!grepl("wavefetch", nams)]
  
  env_layers <- terra::subset(env_layers$layers, nams)
  
  # env_layers <- crop(env_layers, clip_area)
  # 
  # env_layers <- mask(env_layers, mask_area, inverse = T)
  
  vif_list[[i]] <- usdm::vif(env_layers)
}

# Print the results
for (i in to_check) {
  cat(i, ">>", test_grid[i, 1], test_grid[i, 2], "\n")
  print(vif_list[[i]])
  cat("\n")
}


# The problem is, in general, between thetao-* and some other variable
# The high correlation between SST and air temperature was expected.
# However, because patterns change specially close to the poles, we will keep both

# We will now try to resolve the problem between rugosity and slope

# To open the doc:
# rstudioapi::documentOpen("sdm_conf.yml")


# After edited, check again
for (i in to_check) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 2],
                               load_all = T,
                               verbose = FALSE)
  
  nams <- unique(unlist(env_layers$hypothesis))
  nams <- nams[!grepl("wavefetch", nams)]
  
  env_layers <- terra::subset(env_layers$layers, nams)
  
  # env_layers <- crop(env_layers, clip_area)
  # 
  # env_layers <- mask(env_layers, mask_area, inverse = T)
  
  vif_list[[i]] <- usdm::vif(env_layers)
}

# Print the results
for (i in to_check) {
  cat(i, ">>", test_grid[i, 1], test_grid[i, 2], "\n")
  print(vif_list[[i]])
  cat("\n")
}

# Temperature still have a high VIF. This is expected, temperature is usually
# strongly correlated with other variables.
env_layers <- get_envofgroup(group = test_grid[1, 1],
                             depth = test_grid[1, 2],
                             load_all = T,
                             verbose = T)

nams <- unique(unlist(env_layers$hypothesis))
nams <- nams[!grepl("wavefetch", nams)]
nams <- nams[!grepl("tas_mean", nams)] # Remove just for testing

env_layers <- terra::subset(env_layers$layers, nams)

usdm::vif(env_layers)
check_res <- usdm::vifcor(env_layers)
check_res@corMatrix

# However, SST is an essential variable. We then keep this variable and the remaining for this first run.
# We will also keep O2 and PAR, although they are causing some problems with the VIF

# We save the results
names(vif_step_list) <- paste0(test_grid[,1], "_", test_grid[,2])
names(vif_list) <- names(vif_step_list)[min(to_check):max(to_check)]
all_results <- list(
  vif_before = vif_step_list,
  vif_after = vif_list[-which(sapply(vif_list, is.null))]
)

fs::dir_create("data/log")
saveRDS(all_results, file = paste0("data/log/vif_list_", format(Sys.Date(), "%Y%m%d"), ".rds"))
