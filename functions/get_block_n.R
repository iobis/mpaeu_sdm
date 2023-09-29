get_blocks_n <- function(dat, block_type = "spatial_grid") {
  
  k <- dat$blocks$n
  
  n_points <- data.frame(train_pres = rep(NA, k), train_back = NA,
                         test_pres = NA, test_back = NA)
  
  bblock <- list(fold = dat$blocks$folds[[block_type]],
                 pres = dat$training$presence)
  
  for (z in 1:k) {
    n_points$train_pres[z] <- length(bblock$pres[bblock$fold != z & bblock$pres == 1])
    n_points$train_back[z] <- length(bblock$pres[bblock$fold != z & bblock$pres == 0])
    n_points$test_pres[z] <- length(bblock$pres[bblock$fold == z & bblock$pres == 1])
    n_points$test_back[z] <- length(bblock$pres[bblock$fold == z & bblock$pres == 0])
  }
  
  return(n_points)
}

#get_blocks_n(sp_data, "spatial_lat")
