# download meterological data from NASA power API
# downloaded 4 April 2025
# version of NASA power API
# at https://power.larc.nasa.gov/api/pages/
# is v2.6.9

# active cells -----------------------------------------------------------------

# Lambert azimuthal equal area
# see Steinwand et al. (1995)
a.crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# active cells
library(terra)
active.r <- rast("../covariates_spatial/activeAfrica.tif")
same.crs(crs(active.r), a.crs) # TRUE

# project to WGS84
wgs.r <- project(active.r, "epsg:4326")

# bounding box for active cells
ext(wgs.r)
res(wgs.r)

# NASA Power / MERRA 2 raster --------------------------------------------------
# MERRA-2 underlies the historical data available through NASA POWER
# as described here https://power.larc.nasa.gov/docs/methodology/
# Grid structure for MERRA-2 described here:
# https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf

# algorithm for grid
# lat_i = -180 + delta.lat*(i - 1), i = 1, ..., 576
# lon_j = -90  + delta.lon*(j - 1), j = 1, ..., 361
# check algorithm
delta.lon <- 5/8
delta.lat <- 1/2
lons <- seq(-180, 180 - delta.lon, length = 576); unique(diff(lons)) == delta.lon
lats <- seq(-90, 90, length = 361); unique(diff(lats)) == delta.lat
all(c(lons[289], lats[181]) == c(0, 0)) # check origin

# global raster
m2 <- rast(nrows = 361, ncols = 576,
  xmin = -180 - delta.lon/2, xmax = 180 - delta.lon/2,
  ymin = -90 - delta.lat/2, ymax = 90 + delta.lat/2)
all(res(m2) == c(delta.lon, delta.lat)) # check res
all(origin(m2) - c(delta.lon/2, delta.lat/2) == c(0, 0)) # check raster origin
crds(m2)[rowSums(crds(m2) == 0) == 2, ] # centroid origin
# set crs
crs(m2) <- "epsg:4326"

# set extent to active.r projected to WGS84
m2c <- crop(m2, wgs.r)
# assign value of one if MERRA-2 cell includes at least one active cell
m2c.r <- resample(wgs.r, m2c, method = "sum")
m2c.r <- ifel(m2c.r > 0, m2c.r, NA)

# save raster
saveRDS(m2c.r, file = "MERRA2raster.rds")

# coordinates of MERRA-2 cells
m2.crds <- crds(m2c.r)

m2.df <- as.data.frame(m2c.r, xy = TRUE, cells = TRUE)
m2.df$download.id <- 1:nrow(m2.crds)
write.csv(m2.df, file = "MERRA-2-coordinates.csv")

library(nasapower)

# first cell:
# NASA/POWER CERES/MERRA2 Native Resolution Daily Data
# Dates (month/day/year): 06/01/2001 through 08/01/2021
# Location: Latitude  23.5   Longitude 28.75
# Elevation from MERRA-2: Average for 0.5 x 0.625 degree lat/lon region = 281.36 meters
# The value for missing source data that cannot be computed or is outside of the sources availability range: NA
# Parameter(s):
#
#   Parameters:
#   T2M             MERRA-2 Temperature at 2 Meters (C) ;
# RH2M            MERRA-2 Relative Humidity at 2 Meters (%) ;
# PRECTOTCORR     MERRA-2 Precipitation Corrected (mm/day)

pt <- proc.time()

start.i <- 1

for (i in start.i:nrow(m2.crds)) {
  
  try(
    met <- get_power(
      "ag",
      pars = c("T2M", "RH2M", "PRECTOTCORR"),
      temporal_api = "daily",
      lonlat = c(m2.crds[i, 1], m2.crds[i, 2]),
      dates = c("2001-06-01", "2021-08-01") # local solar time (LST)
    )
  )
  
  Sys.sleep(runif(1, 1.25, 4.75))
  flush.console()
  
  iter <- 1
  while(inherits(met, "try-error") || all(as.matrix(met[1, c(1, 2)]) == as.matrix(met.old[1, c(1, 2)]))) {
 
    cat(iter, "\n")
    
    Sys.sleep(runif(1, 5.5, 10))
    flush.console()
 
    try(
      met <- get_power(
        "ag",
        pars = c("T2M", "RH2M", "PRECTOTCORR"),
        temporal_api = "daily",
        lonlat = c(m2.crds[i, 1], m2.crds[i, 2]),
        dates = c("2001-06-01", "2021-08-01") # local solar time (LST)
      )
    )
    if (iter > 10) stop("too many tries")
    iter <- iter + 1
  }
  
  saveRDS(met, file = paste0("met_", i,".rds"))
  met.old <- met

  cat("completed ", i, " of ", nrow(m2.crds), "\n")
  cat((proc.time()['elapsed'] - pt['elapsed'])/60, " min \n")
    
}

