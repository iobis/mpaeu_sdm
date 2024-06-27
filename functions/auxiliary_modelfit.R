# Internal functions - to document

# Retrieve any object from the results
retrieve <- function(species, what, acro = NULL, results_folder = "results/",
                     load = TRUE, return_all = FALSE) {
  
  if (return_all) {
    load <- FALSE
  }
  
  all_f <- list.files(results_folder)
  
  all_f <- all_f[grepl(species, all_f)]
  
  if (length(all_f) < 1) {
    stop("No folder found for that species.")
  }
  
  all_files <- list.files(paste0(results_folder, "/", all_f), full.names = T)
  
  if (!is.null(acro)) {
    all_files <- all_files[grepl(acro, all_files)]
    if (length(all_files) < 1) {
      stop("No folder found for this acronym.")
    }
  } else {
    folders_info <- fs::dir_info(paste0(results_folder, "/", all_f))
    all_files <- unlist(folders_info[which.max(folders_info$modification_time), "path"])
  }
  
  all_files <- list.files(all_files, recursive = T, full.names = T)
  
  if (length(all_files) < 1) {
    stop("Directory is empty.")
  }
  
  sel_file <- all_files[grepl(what, all_files)]
  
  if (length(sel_file) > 1) {
    if (!return_all) {
      if (interactive()) {
        print(sel_file)
        chose <- as.numeric(readline("Which file you want?  "))
        if (!chose %in% 1:length(sel_file)) {
          cat("File must be one of the available. Supply a number!\n")
          chose <- as.numeric(readline("Which file you want?  "))
        }
      } else {
        cat("Multiple files, returning the first.\n")
        chose <- 1
      }
      sel_file <- sel_file[chose]
    }
  }
  
  
  if (load) {
    result <- switch (tools::file_ext(sel_file[1]),
                      parquet = arrow::read_parquet(sel_file),
                      tif = terra::rast(sel_file),
                      shp = terra::vect(sel_file),
                      csv = read.csv(sel_file),
                      rds = readRDS(sel_file),
                      {cat("Impossible to load file with extension", tools::file_ext(sel_file), "\n");NULL}
    )
  } else {
    result <- sel_file
  }
  
}
  