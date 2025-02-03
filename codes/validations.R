############################## MPA Europe project ##############################
########### WG3 - Species and biogenic habitat distributions (UNESCO) ##########
# September of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################ C++ codes validation ##############################

# Classification
library(terra)

Rcpp::sourceCpp("functions/classify_raster.cpp")

test_r <- rast(nrow = 5, ncol = 5, nlyr = 4)
values(test_r) <- rep(seq(30, 79), 2)

plot(test_r)
minmax(test_r)

mins <- minmax(test_r)[1,]

ths <- mins + 5

tictoc::tic()
classified <- app(test_r, applyThreshold, thresholds = ths)
tictoc::toc()

plot(classified)

classified[classified == 0] <- NA

minmax(test_r)
minmax(classified)
cat("Expected:", ths)
cat("Get:", minmax(classified)[1,])

par(mfrow = c(2,2))
plot(test_r[[1]])
plot(classified[[1]], main = ths[1])
plot(test_r[[2]])
plot(classified[[2]], main = ths[2])


### END