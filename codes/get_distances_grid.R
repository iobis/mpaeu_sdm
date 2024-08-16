############################## MPA Europe project ##############################
########### WP3 - Species and biogenic habitat distributions (UNESCO) ##########
# March of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
#
########################### Generate distance grids ############################

library(obissdm)
library(terra)
library(biooracler)

base <- download_layers("thetao_baseline_2000_2019_depthsurf",
                        "thetao_mean", list(
                            time = c('2010-01-01T00:00:00Z', '2010-01-01T00:00:00Z'),
                            latitude = c(-89.975, 89.975),
                            longitude = c(-179.975, 179.975)
                        ))

base
base[!is.na(base)] <- 1
plot(base)

outqc_get_distances(base = base,
                    agg_res = 0.8,
                    outfolder = "data/distances",
                    mc_cores = 4)

outqc_dist_tords("data/distances", parallel_cores = 8, remove_distances = T)

###END