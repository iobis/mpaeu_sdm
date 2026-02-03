library(duckdb)
library(DBI)
library(h3jsr)
library(terra)
library(spatstat)
library(sf)

con <- dbConnect(duckdb())

send <- function(...) {dbSendQuery(con, ...)}
get <- function(...) {dbGetQuery(con, ...)}

send("install h3 from community; load h3;")
send("install icu; load icu;")

effort_total <- get(
    "select 
        h3_latlng_to_cell_string(decimalLatitude, decimalLongitude, 6) as h3_6,
        count(distinct(
            concat(datepart('year', to_timestamp(date_mid / 1000)), '-',
               lpad(cast(datepart('month', to_timestamp(date_mid / 1000)) as varchar), 2, '0'))
        )) as unique_months
    from read_parquet('data/raw/obis_20240625.parquet')
    where date_mid is not null
    group by h3_6
    "
)

effort_total$normalized_effort <- (effort_total$unique_months - min(effort_total$unique_months)) /
    (max(effort_total$unique_months) - min(effort_total$unique_months))

effort_total_xy_s <- sf::st_as_sf(effort_total_xy, coords = c("x", "y"), crs = "EPSG:4326")

win <- as.owin(st_bbox(effort_total_xy_s))  # Define window extent
points_ppp <- ppp(x = st_coordinates(effort_total_xy_s)[,1],
                  y = st_coordinates(effort_total_xy_s)[,2],
                  window = win, marks = effort_total_xy_s$effort)

density_map <- Smooth(points_ppp, sigma = 5)

plot(rast(density_map))
samp <- sample(seq_len(nrow(effort_total_xy)), 10000)
points(effort_total_xy[samp,1:2], cex = effort_total_xy$effort[samp]*10)

density_map_in <- rast(density_map)
crs(density_map_in) <- crs("EPSG:4326")
density_map_in <- project(density_map_in, rast(res = 0.05, crs = "EPSG:4326"))

base_raster <- rast("data/env/current/thetao_baseline_depthsurf_mean.tif")
density_map_in <- mask(density_map_in, base_raster)
plot(density_map_in)
writeRaster(density_map_in, "data/virtual_species/sampling_effort.tif")

#### Other ideas
# Perform kernel density estimation (adjust bandwidth for smoothing)
# density_map <- density(points_ppp, sigma = 5)  # Sigma controls smoothing
# Generate a global map of effort
# effort_total_cord <- cell_to_point(effort_total$h3_6)
# effort_total_xy <- as.data.frame(sf::st_coordinates(effort_total_cord))
# colnames(effort_total_xy) <- c("x", "y")
# effort_total_xy$effort <- effort_total$normalized_effort

# v <- variogram(effort ~ 1, ~x+y, data = effort_total_xy[sample(seq_len(nrow(effort_total_xy)), 100000),])
# mv <- fit.variogram(v, vgm(1, "Sph", 300, 1))


# library(mgcv)

# smooth_gam <- bam(effort ~ s(x, y, bs = 'gp', k = 100, m = 2), 
#                   data = effort_total_xy)

# base_r <- rast(res = 0.5)
# base_r[] <- 1
# base_r_cd <- as.data.frame(base_r, xy = T)

# base_r_res <- predict(smooth_gam, base_r_cd)

# base_r[cellFromXY(base_r, base_r_cd[,1:2])] <- as.vector(base_r_res)

# plot(base_r)
# samp <- sample(seq_len(nrow(effort_total_xy)), 10000)
# points(effort_total_xy[samp,1:2], cex = effort_total_xy$effort[samp]*10)

# library(sf)



# r_smooth <- terra::density(vect(effort_total_xy_s), field = "effort", radius = 5)



# idw_v <- interpIDW(rast(res = 0.5), vect(effort_total_xy_s), "effort", radius = 20, maxPoints = 30)
# plot(idw_v)

# library(fields) 
# samp <- sample(seq_len(nrow(effort_total_xy)), 10000)
# tps <- Tps(as.matrix(effort_total_xy[samp,1:2]), effort_total_xy$effort[samp], lon.lat = T)
# p <- rast(r)
# p <- interpolate(p, tps)
# p <- mask(p, r)
# plot(p)

# mg <- gstat(id = "effort", formula = effort~1, locations = ~x+y, data=effort_total_xy, 
#             nmax=15)
# z <- interpolate(rast(res = 0.5), mg)
# plot(z)
