ifok <- function(result, ...){
  error_message <- NA
  if (class(result)[1] != "try-error") {
    eval(...)
  } else {
    error_message <- attr(result, "condition")$message
  }
  #return(invisible(NULL))
  return(error_message)
}

save_model_stuff <- function(model, 
                             env_layers,
                             outfile) {
  
  # Save model
  saveRDS(model, file = paste0(outfile, "_model.rds"))
  
  # Save metrics
  write.csv(model$cv_metrics, paste0(outfile, "_cvmetrics.csv"), row.names = F)
  full_metrics <- data.frame(
    metrics = names(model$full_metrics),
    values = model$full_metrics
  )
  write.csv(full_metrics, paste0(outfile, "_fullmetrics.csv"), row.names = F)
  
  if (!is.null(model$eval_metrics)) {
    write.csv(model$eval_metrics, paste0(outfile, "_evalmetrics.csv"), row.names = F)
  }
  
  # Save response curves
  respc_data <- resp_curves(model, env_layers)
  write.csv(respc_data, paste0(outfile, "_respcurves.csv"), row.names = F)
  respc_plot <- plot(respc_data)
  ggsave(paste0(outfile, "_respcurves.jpg"), respc_plot, quality = 100,
         width = 6, height = 4, units = "in")
  
  return(invisible(NULL))
}


log_start <- function(outfile, ...) {
  
  if (length(list(...)) > 0) {
    env <- environment()
    list2env(list(...), envir = env)
  }
  
  fs::dir_create(dirname(outfile))
  
  # Save log file
  to_write <- c(
    "# LOG FILE FOR SDM MODELING - MPA EUROPE #",
    "##########################################",
    paste0(outacro, ":"),
    paste("  date_first_run:", format(Sys.Date(), "%Y-%m-%d")),
    #paste("  algos:", paste(algos, collapse = ";")),
    paste("  minimum_pts:", min_n),
    paste("  tried_eval_iflown:", iflow_tryeval),
    paste("  behavior_with_higheval:", ifhigheval),
    paste("  help_perc:", help_perc),
    paste("  standard_habitat:", std_habitat),
    paste("  try_other_blocks:", try_other_blocks),
    paste("  min_size_block:", min_size_block),
    paste("  max_size_block:", max_size_block),
    paste("  clean_sac:", min_size_block),
    paste("  max_dist_sac:", max_sac),
    paste("  hypothesis_tested:", paste(test_multi_hypo, collapse = ";")),
    paste("  future_scenarios:", paste(fut_scen, collapse = ";")),
    paste("  future_periods:", paste(fut_periods, collapse = ";")),
    paste("  clip_area:", ifelse(is.null(clip_area), NA,
                                 paste(as.vector(clip_area), collapse = ";"))),
    paste("  limited_depth:", limit_by_depth),
    paste("  variables_scaled:", scale_vars),
    paste("  species:")
  )
  
  writeLines(to_write, con = outfile)
  
  outfile <- paste0(dirname(outfile), "/_logtemp")
  fs::dir_create(outfile)
  
  return(outfile)
  #return(invisible(NULL))
}

log_addsp <- function(log_file, sp, redo = FALSE) {
  
  #prev <- readLines(log_file, warn = F)
  outfile <- paste0(log_file, "/", sp, "_temp.yml")
  
  if (!file.exists(outfile) | redo) {
    to_write <- c(
      #prev,
      # Level 3 (4 spaces)
      paste0("    key", sp, ":")
    )
    
    writeLines(to_write, con = outfile)
  }
  
  #return(invisible(NULL))
  return(outfile)
}

log_addmodel <- function(log_file, model_name, model, sp_data, sp_max_depth,
                         tune_method, block_size, clean_up_mode = TRUE) {
  
  prev <- readLines(log_file, warn = F)
  
  if (clean_up_mode) {
    if (any(grepl(model_name, prev))) {
      prev <- prev[-grep(model_name, prev):-(grep(model_name, prev)+9)]
    }
  }
  
  to_write <- c(
    prev,
    # Level 4 (6 spaces)
    paste0("      ", model_name, ":"),
    # Level 5 (8 spaces)
    paste("        variables_used:", paste0(colnames(sp_data$training)[-1], collapse = ";")),
    paste("        depth:", sp_max_depth),
    paste("        n_presences:", sum(sp_data$training$presence)),
    paste("        n_quadrature:", length(sp_data$training$presence[sp_data$training$presence == 0])),
    paste("        n_evaluation:", ifelse(is.null(sp_data$eval_data), 0, sum(sp_data$eval_data$presence))),
    paste("        tune_method:", tune_method),
    paste("        cv_method:", paste0(model$cv_method, collapse = ";")),
    paste("        final_block_size:", block_size),
    paste("        timings:", paste0(paste(names(model$timings), round(model$timings, 2), sep = "-"), collapse = ";"))
  )
  
  writeLines(to_write, con = log_file)
  
  return(invisible(NULL))
}

log_close <- function(log_file) {
  
  basedir <- dirname(log_file)
  
  temp_folder <- paste0(basedir, "/_logtemp")
  
  all_files <- list.files(temp_folder, full.names = T)
  
  init_file <- readLines(log_file)
  
  for (i in 1:length(all_files)) {
    temp_file <- readLines(all_files[i])
    
    if (length(temp_file) > 2) {
      init_file <- c(init_file, temp_file)
    }
  }
  
  writeLines(init_file, con = log_file)
  
  cli::cli_alert_success("Log file {.file {log_file}} successfuly closed.")
  
  fs::dir_delete(temp_folder)
  
  return(invisible(NULL))
}
