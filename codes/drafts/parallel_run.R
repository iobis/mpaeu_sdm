parallel_modeling <- function(
    species,
    eco_info,
    sp_groups,
    log_file,
    #flux_ctrl_name,
    outfold,
    outacro,
    skip_done,
    algos,
    run_mode,
    min_n,
    iflow_tryeval,
    ifhigheval,
    help_perc,
    clean_sac,
    max_sac,
    std_habitat,
    try_other_blocks,
    min_size_block,
    max_size_block,
    test_multi_hypo,
    fut_scen,
    fut_periods,
    clip_area,
    mask_area_names,
    limit_by_depth,
    depth_buffer,
    depth_layer,
    scale_vars
){
  
  set.seed(2023)
  
  #flux_ctrl <- storr_rds(flux_ctrl_name)
  
  if (!is.null(mask_area_names)) {
    mregions <- mregions::mr_shp("MarineRegions:iho")
    mask_area <- mregions[mregions$name %in% mask_area_names,]
  }
  
  depth_layer <- rast(depth_layer)
  
  all_sp <- sp_groups
  
  run <- run_models(species, all_params = mget(ls()))
  
  return(invisible(NULL))
}


run_parallel <- function(cluster, species_list, ...) {
  
  if (!file.exists("data/log/sdm_models_log.yml")) {
    log_file <- log_start("data/log/sdm_models_log.yml", ...)
  } else {
    log_file <- "data/log/_logtemp"
  }
  # 
  # cl <- makeCluster(ncores)
  result <- parLapply(cluster, species_list, parallel_modeling,
                      log_file = log_file, ...)
  # stopCluster(cl)
  # result <- lapply(species_list, parallel_modeling,
  #                           log_file = log_file, ...)
  
  return(invisible(NULL))
  
}
