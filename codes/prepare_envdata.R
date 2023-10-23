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

mregions <- mregions::mr_shp("MarineRegions:iho")
mask_area <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

groups <- names(get_conf(what = "groups")$groups)
depth_env <- c("depthsurf", "depthmean", "depthmax")

test_grid <- expand.grid(groups, depth_env, stringsAsFactors = F)


# Create a list to hold the results of vif step
vif_step_list <- list()

for (i in 1:nrow(test_grid)) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 2],
                               verbose = TRUE)
  
  env_layers <- crop(env_layers, clip_area)
  
  env_layers <- mask(env_layers, mask_area, inverse = T)
  
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
                               verbose = FALSE)
  
  env_layers <- crop(env_layers, clip_area)
  
  env_layers <- mask(env_layers, mask_area, inverse = T)
  
  vif_list[[i]] <- usdm::vif(env_layers)
}

# Print the results
for (i in to_check) {
  cat(i, ">>", test_grid[i, 1], test_grid[i, 2], "\n")
  print(vif_list[[i]])
  cat("\n")
}


# The problem is, in general, between thetao-* and some other variable
# We will try to remove those suggested by VIF step or keep if temperature
# To open the doc:
# rstudioapi::documentOpen("sdm_conf.yml")


# After edited, check again
for (i in to_check) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 2],
                               verbose = FALSE)
  
  env_layers <- crop(env_layers, clip_area)
  
  env_layers <- mask(env_layers, mask_area, inverse = T)
  
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
env_layers <- get_envofgroup(group = test_grid[2, 1],
                             depth = test_grid[2, 2],
                             verbose = T)

env_layers <- crop(env_layers, clip_area)

env_layers <- mask(env_layers, mask_area, inverse = T)

check_res <- usdm::vifcor(env_layers)
check_res@corMatrix

# However, SST is an essential variable. We then keep this variable and the remaining for this first run.

# We save the results
names(vif_step_list) <- paste0(test_grid[,1], "_", test_grid[,2])
names(vif_list) <- names(vif_step_list)[1:11]
all_results <- list(
  vif_before = vif_step_list,
  vif_after = vif_list[-which(sapply(vif_list, is.null))]
)

fs::dir_create("data/log")
saveRDS(all_results, file = paste0("data/log/vif_list_", format(Sys.Date(), "%Y%m%d"), ".rds"))
