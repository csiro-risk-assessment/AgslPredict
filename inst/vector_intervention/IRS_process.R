# This script process IRS data downloaded from # https://malariaatlas.org/
# about Africa IRS Coverage corresponding to the 
# proportion of households covered with Indoor Residual Spraying during a 
# defined year
# The script reads the TIFF files 2024_GBD2023_Africa_IRS_XXXX.tif where 
# XXXX goes from 2000 to 2021 and projects the raster on the active grid then
# saves the raster into TIFF files

library(terra)

# reading the active grid used to project the population
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
active.r <- rast("../covariates_spatial/activeAfrica.tif")

# please set the variable IRS.dir to the path name of the directory
# where the files downloaded from https://malariaatlas.org/ are located
# you need the TIFF files 2024_GBD2023_Africa_IRS_XXXX.tif where 
# XXXX goes from 2000 to 2021
IRS.dir <- " "

# project to WGS84
wgs.r <- project(active.r, "epsg:4326")

years = seq(2000,2021)

for (year in years)
{
  filename = paste0(IRS.dir,"2024_GBD2023_Africa_IRS_", year, ".tif")
  IRS.r <- rast(filename, win = ext(wgs.r))
  pIRS <- project(IRS.r, crs(active.r))
  IRS.r <- resample(pIRS, active.r, method = 'average') 
  writeRaster(IRS.r, filename = paste0("./IRS_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff")
}

