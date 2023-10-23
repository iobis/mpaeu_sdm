library(h3)
h3_add <- geo_to_h3(species_points[,c("decimalLongitude", "decimalLatitude")], res = 3)

teste <- h3_to_geo_boundary_sf(h3_add)
plot(teste)
teste$id <- 1:nrow(teste)

mapview::mapview(teste[,"id"])

teste_b <- teste[c(50,30,48,36,46,49,16,18),]

dist_vals <- rep(NA, nrow(teste_b))

target_pt <- teste_b[8,"h3_index"] %>%
  sf::st_drop_geometry() %>%
  unlist

for (i in nrow(teste_b):1) {
  
  start_pt <- teste_b[i,"h3_index"] %>%
    sf::st_drop_geometry() %>%
    unlist()
  
  route <- start_pt
  
  while (start_pt != target_pt) {
    nearby <- k_ring(start_pt, radius = 1)
    closest <- which.min(h3_distance(target_pt, nearby))[1]
    start_pt <- nearby[closest]
    route <- c(route, start_pt)
    print("1st.\n")
  }
  
  
}

dist_matrix <- matrix(nrow = nrow(teste_b), ncol = 1)

not_valid <- c("837b69fffffffff", "837b6bfffffffff", "837b6dfffffffff")
for (k in 1:nrow(teste_b)) {
  
  start_pt <- teste_b[k,"h3_index"] %>%
    sf::st_drop_geometry() %>%
    unlist()
  #start_pt <- "8396b4fffffffff"
  # 
  # route <- start_pt
  # 
  if (!start_pt %in% not_valid) {
    rad <- 1
    nearby <- c()
    while (!target_pt %in% nearby) {
      nearby <- k_ring(start_pt, radius = rad)
      if (!target_pt %in% nearby) {
        rad <- rad+1
      }
      print("1st.\n")
    }
    
    valid <- nearby[!nearby %in% not_valid]
    
    original_start <- route <- start_pt
    
    # sf_p <- h3_to_geo_boundary_sf(valid)
    # sf_p$route <- 0
    # sf_p$route[sf_p$h3_index == original_start] <- 1
    # plot(sf_p[,"route"])
    
    while (start_pt != target_pt) {
      nearby_c <- k_ring(start_pt, radius = 1)
      nearby_c <- nearby_c[!nearby_c %in% original_start]
      nearby_c <- nearby_c[nearby_c %in% valid]
      closest_dists <- h3_distance(target_pt, nearby_c)
      closest <- which(closest_dists == min(closest_dists))
      ### Check
      # sf_p$route <- ifelse(sf_p$h3_index %in% nearby_c[closest], 2, sf_p$route)
      # plot(sf_p[,"route"])
      ###
      min_dist <- c()
      i <- 1
      while (i < (length(closest)+1)) {
        nearby_t <- k_ring(nearby_c[closest[i]], radius = 1)
        nearby_t <- nearby_t[!nearby_t %in% c(original_start, nearby_c[closest[i]])]
        nearby_t <- nearby_t[nearby_t %in% valid]
        closest_dists_t <- h3_distance(target_pt, nearby_t)
        min_dist <- c(min_dist, min(closest_dists_t))
        i <- i+1
      }
      start_pt <- nearby_c[closest[which.min(min_dist)]]
      # sf_p$route[sf_p$h3_index == start_pt] <- 3
      # plot(sf_p[,"route"])
      #start_pt <- nearby_c[closest]
      route <- c(route, start_pt)
      print("1st.\n")
      
      # sf_p <- h3_to_geo_boundary_sf(valid)
      # sf_p$route <- ifelse(sf_p$h3_index %in% route, 1, 0)
      # plot(sf_p)
      
      #Sys.sleep(2)
    }
    
    dist_matrix[(nrow(dist_matrix) - k)+1,] <- length(route)
    # sf_p <- h3_to_geo_boundary_sf(valid)
    # sf_p$route <- ifelse(sf_p$h3_index %in% route, 1, 0)
    #plot(sf_p[,"route"], main = "final")
    #Sys.sleep(4)
  }
}



