# Calculate the ideal number of quadrature points

n.quad <- c(seq(5000, 30000, by = 5000), seq(30000, 50000, by = 10000),
            seq(50000, 500000, by = 50000))

for (key in 101:104) {
  lik.list <- list()
  
  sp_key <- key
  nsamp <- "low"
  
  train_data <- read_parquet(paste0("data/virtual_species/key=", sp_key, "/occurrences_", nsamp, "_rep1.parquet"))
  
  for (j in 1:3) {
    
    lik <- rep(NA, length(n.quad))
    
    for (i in 1:length(n.quad)) {
      
      sp_data <- mp_prepare_data(train_data, species_id = paste0("species", sp_key),
                                 env_layers = curr,
                                 quad_number = n.quad[i]
                                 #, pred_quad = pred_quad_pts
      )
      
      p.wt <- rep(1.e-8, n.quad[i])
      p.wt[sp_data$training$presence == 0] <- dwpr_val/n.quad[i]
      
      z <- sp_data$training$presence/p.wt
      
      dwpr <- glm(z ~ thetao_mean + so_mean + po4_mean + phyc_mean + bathymetry_mean,
                  family = poisson(), weights = p.wt, data = sp_data$training)
      
      mu <- dwpr$fitted
      
      lik[i] <- sum(p.wt*(z*log(mu) - mu))
    }
    
    lik.list[[j]] <- lik
    
    cat(j, "|")
  }
  
  plot(y = lik.list[[1]], x = n.quad, type = "o", log = "x",
       ylim = range(unlist(lik.list)))
  lapply(2:length(lik.list), function(x){
    lines(y = lik.list[[x]], x = n.quad, type = "o", col = sample(colors(), 1))
  })
}