### See the ideal block size for spatial cross-validation
library(terra)
set.seed(2023)
sf::sf_use_s2(FALSE)

clip_area <- ext(-41, 47, 20, 89)

mregions <- mregions::mr_shp("MarineRegions:iho")
mask_area <- mregions[mregions$name %in% c("Red Sea", "Gulf of Aqaba", "Gulf of Suez"),]

groups <- names(get_conf(what = "groups")$groups)
depth_env <- c("depthsurf", "depthmean", "depthmax")

test_grid <- expand.grid(groups, depth_env, stringsAsFactors = F)
test_grid$size <- NA

for (i in 1:nrow(test_grid)) {
  env_layers <- get_envofgroup(group = test_grid[i, 1],
                               depth = test_grid[i, 1],
                               verbose = FALSE)
  
  env_layers <- crop(env_layers, clip_area)
  
  env_layers <- mask(env_layers, mask_area, inverse = T)
  
  autocor <- blockCV::cv_spatial_autocor(env_layers,
                                         plot = F)
  test_grid$size[i] <- autocor$range/1000
}

test_grid$size <- round(test_grid$size, 2)

conf_file <- readLines("sdm_conf.yml")
conf_file <- c(
  conf_file, "# Cross-validation block sizes",
  "# Block size assessed using the package blockCV",
  "# See the code get_block_size.R for more info",
  "blocksizes:",
  "  depthsurf:",
  paste0("    ", test_grid[test_grid[,2] == "depthsurf",1], ": ",
         test_grid[test_grid[,2] == "depthsurf",3]),
  "  depthmean:",
  paste0("    ", test_grid[test_grid[,2] == "depthmean",1], ": ",
         test_grid[test_grid[,2] == "depthmean",3]),
  "  depthmax:",
  paste0("    ", test_grid[test_grid[,2] == "depthmax",1], ": ",
         test_grid[test_grid[,2] == "depthmax",3])
)
writeLines(conf_file, con = "sdm_conf.yml")
