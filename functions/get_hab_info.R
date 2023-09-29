get_hab_info <- function(species,
                         databases = c("sealife", "fishbase", "worms"),
                         #biotic_file = NULL,
                         try_higher = TRUE,
                         search_remarks = TRUE) {
  
  if ("biotic" %in% databases) {
    if (is.null(biotic_file)) {
      stop("For BIOTIC, a file containing traits is needed.")
    }
  }
  
  worms <- fishbase <- sealife <- "NOT_FOUND"
  #biotic <- NULL
  
  species_info <- worrms::wm_record(species)
  species_info <- species_info[species_info$status == "accepted",]
  species_info <- species_info[1,]
  
  species_info <- species_info[,c("AphiaID", "scientificname", "genus", "family", "order", "class")]
  
  if ("sealife" %in% databases) {
    sealife <- .info_from_sealife(species_info$scientificname)
  }
  if ("fishbase" %in% databases) {
    fishbase <- .info_from_sealife(species_info$scientificname, server = "fishbase")
  }
  if ("worms" %in% databases) {
    worms <- .info_from_worms(species)
  }
  
  if (sealife != "NOT_FOUND") {
    result <- sealife
  } else if (fishbase != "NOT_FOUND") {
    result <- fishbase
  } else if (worms != "NOT_FOUND") {
    result <- worms
  } else if (try_higher) {
    result <- .info_from_worms(species, try_parent = TRUE)
  } else {
    result <- "NOT_FOUND"
  }
  
  return(result)
  
}


.info_from_sealife <- function(species_name, server = "sealife") {
  
  info <- rfishbase::ecology(species_name, server = server)[1,]
  colnames(info) <- tolower(colnames(info))
  
  if (nrow(info) < 1) {
    mode_life <- "NOT_FOUND"
  } else {
    test_info <- function(colname) {
      if (colname %in% colnames(info)) {
        ifelse(info[1, colname] == 1, TRUE, FALSE)[,1]
      } else {
        NULL
      }
    }
    
    is_benthic <- test_info("benthic")
    
    is_demersal <- test_info("demersal")
    
    is_pelagic <- test_info("pelagic")
    
    mode_life <- ifelse(isTRUE(is_benthic), "benthic",
                        ifelse(isTRUE(is_demersal), "demersal",
                               ifelse(isTRUE(is_pelagic), "pelagic", "NOT_FOUND")))
    
    if (mode_life == "pelagic") {
      is_epipelagic <- test_info("epipelagic")
      
      is_mesopelagic <- test_info("mesopelagic")
      
      is_bathipelagic <- test_info("bathypelagic")
      is_abyssopelagic <- test_info("abyssopelagic")
      is_hadopelagic <- test_info("hadopelagic")
      
      mode_life <- ifelse(isTRUE(is_epipelagic), "pelagic_surface",
                          ifelse(isTRUE(is_mesopelagic), "pelagic_mean",
                                 ifelse(isTRUE(is_bathipelagic) | isTRUE(is_abyssopelagic) | isTRUE(is_hadopelagic),
                                        "pelagic_bottom", "pelagic_surface")))
    }
    
    if (mode_life == "NOT_FOUND") {
      sp_info <- rfishbase::species(species_name, server = server)
      colnames(sp_info) <- tolower(colnames(sp_info))
      
      if ("demerspelag" %in% colnames(sp_info)) {
        is_benthic <- ifelse(sp_info$demerspelag == "benthic", TRUE, FALSE)
        is_demersal <- ifelse(sp_info$demerspelag == "demersal", TRUE, FALSE)
        is_pelagic <- ifelse(sp_info$demerspelag == "pelagic", TRUE, FALSE)
        
        mode_life <- ifelse(isTRUE(is_benthic), "benthic",
                            ifelse(isTRUE(is_demersal), "demersal",
                                   ifelse(isTRUE(is_pelagic), "pelagic", "NOT_FOUND")))
      }
      
      if (mode_life == "NOT_FOUND") {
        if ("comments" %in% colnames(sp_info)) {
          com_info <- tolower(sp_info$comments)
          
          is_benthic <- ifelse(grepl("benthic", com_info), TRUE, FALSE)
          is_demersal <- ifelse(grepl("demersal", com_info), TRUE, FALSE)
          is_pelagic <- ifelse(grepl("pelagic", com_info), TRUE, FALSE)
          
          mode_life <- ifelse(isTRUE(is_benthic), "benthic",
                              ifelse(isTRUE(is_demersal), "demersal",
                                     ifelse(isTRUE(is_pelagic), "pelagic", "NOT_FOUND")))
        }
      }
    }
  }
  
  return(mode_life)
  
}


.info_from_worms <- function(species_code, try_parent = FALSE) {
  
  info <- worrms::wm_attr_data(species_code, include_inherited = try_parent)
  
  if ("Functional group" %in% info$measurementType) {
    info_life <- info$measurementValue[info$measurementType == "Functional group"]
    
    is_benthic <- any(grepl("benth", tolower(info_life)))
    
    is_demersal <- any(grepl("demers", tolower(info_life)))
    
    is_pelagic <- any(grepl("pelag", tolower(info_life)))
    
    mode_life <- ifelse(isTRUE(is_benthic), "benthic",
                        ifelse(isTRUE(is_demersal), "demersal",
                               ifelse(isTRUE(is_pelagic), "pelagic", "NOT_FOUND")))
    
    if (mode_life == "pelagic") {
      is_epipelagic <- any(grepl("epipelagic", tolower(info_life)))
      
      is_mesopelagic <- any(grepl("mesopelagic", tolower(info_life)))
      
      is_bathipelagic <- any(grepl("bathypelagic", tolower(info_life)))
      is_abyssopelagic <- any(grepl("abyssopelagic", tolower(info_life)))
      is_hadopelagic <- any(grepl("hadopelagic", tolower(info_life)))
      
      mode_life <- ifelse(isTRUE(is_epipelagic), "pelagic_surface",
                          ifelse(isTRUE(is_mesopelagic), "pelagic_mean",
                                 ifelse(isTRUE(is_bathipelagic) | isTRUE(is_abyssopelagic) | isTRUE(is_hadopelagic),
                                        "pelagic_bottom", "pelagic_surface")))
    }
    
  } else {
    mode_life <- "NOT_FOUND"
  }
  
  return(mode_life)
  
}

# .info_from_biotic <- function(species_name, file) {
#   
#   info <- read.csv(file)
#   
#   if (nrow(info) < 1) {
#     mode_life <- "NOT_FOUND"
#   } else {
#     test_info <- function(colname) {
#       if (colname %in% colnames(info)) {
#         ifelse(info[1, colname] == 1, TRUE, FALSE)[,1]
#       } else {
#         NULL
#       }
#     }
#     
#     is_benthic <- test_info("benthic")
#     
#     is_demersal <- test_info("demersal")
#     
#     is_pelagic <- test_info("pelagic")
#     
#     mode_life <- ifelse(isTRUE(is_benthic), "benthic",
#                         ifelse(isTRUE(is_demersal), "demersal",
#                                ifelse(isTRUE(is_pelagic), "pelagic", "NOT_FOUND")))
#     
#     if (mode_life == "pelagic") {
#       is_epipelagic <- test_info("epipelagic")
#       
#       is_mesopelagic <- test_info("mesopelagic")
#       
#       is_bathipelagic <- test_info("bathypelagic")
#       is_abyssopelagic <- test_info("abyssopelagic")
#       is_hadopelagic <- test_info("hadopelagic")
#       
#       mode_life <- ifelse(isTRUE(is_epipelagic), "pelagic_surface",
#                           ifelse(isTRUE(is_mesopelagic), "pelagic_mean",
#                                  ifelse(isTRUE(is_bathipelagic) | isTRUE(is_abyssopelagic) | isTRUE(is_hadopelagic),
#                                         "pelagic_bottom", "pelagic_surface")))
#     }
#   }
#   
#   return(mode_life)
#   
# }