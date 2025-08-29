# This script process the human population data downloaded from
# https://landscan.ornl.gov/ to obtain TIFF files giving the total human
# population and the human density on a 5km by 5km cell, projected on the
# active grid.

library(terra)

# reading the active grid used to project the population
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
active.r <- rast("../covariates_spatial/activeAfrica.tif")

# please set the variable preprocessing.dir to the pathname of the directory
# where the files downloaded from https://landscan.ornl.gov/ are located
# you need the TIFF files landscan-global-XXXX.tif where XXXX goes from 2000 to 2021
preprocessing.dir <- " "
tif_files <- list.files(preprocessing.dir, pattern = "\\.tif$", full.names = TRUE)


# project to WGS84
wgs.r <- project(active.r, "epsg:4326")
print(ext(active.r))
active.mask <- active.r != 0

for (f in tif_files) {
  dem.r <- rast(f, win = ext(wgs.r))
  year <- substring(f, nchar(f)-7, nchar(f)-4)
  print(year)
  # project raster
  dem.proj <- project(dem.r, crs(active.r))
  # resample
  # give the number of human for 25 square km
  dem.resamp.sum <- resample(dem.proj, active.r, method = "sum")
  dem.sum.r <- mask(dem.resamp.sum, active.mask, maskvalues = FALSE)
  # give the density of human by km
  # raster of cell size
  dem.proj.size <- cellSize(dem.proj)
  # we ignore cell with NA pop when calculating the area used in density
  dem.proj.size[is.na(dem.proj)] <- 0
  dem.size <- resample(dem.proj.size, active.r, method = "sum")
  # cell size in km
  dem.size <- dem.size/1e06
  # calculate density
  dem.avg.r <- dem.sum.r/ dem.size
  writeRaster(dem.sum.r, filename = paste0("./human_pop_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff", overwrite = TRUE)
  writeRaster(dem.avg.r, filename = paste0("./human_density_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff", overwrite = TRUE)
}
