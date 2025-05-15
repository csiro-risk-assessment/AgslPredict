# Processing script for data from Malaria Atlas Project
# 202406_Africa_Insecticide_Treated_Net_Use available from Malaria Atlas Project
# URL: https://malariaatlas.org/
# Collection of tiff files for years 2000-2022
# Download date: 31 March 2025

# read in active cells
library(terra)

# Lambert azimuthal equal area
# see Steinwand et al. (1995)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# check active cells for correct extent info
all(read.csv("../PMB_data_2024023/data/processed/Covariates/active.csv",
             header = FALSE, nrow = 1) == c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# read in spatial scope of RA in projected CRS
active <- read.csv("../PMB_data_2024023/data/processed/Covariates/active.csv",
                   skip = 1, header = FALSE)
all.equal(dim(active), c(1280, 1520)) # expecting same number of rows and columns for all covariate data
active <- as.matrix(active)
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

# project to WGS84
wgs.r <- project(active.r, "epsg:4326")

# set data directory for unzipped TIFF files downloaded from Malaria Atlas Project, e.g.:
ITN.dir <- "../2024_GBD2023_Africa_ITN_2000/"

years <- 2000:2022

for (year in years)
{
  filename = paste0(ITN.dir,"2024_GBD2023_Africa_ITN_", year, ".tif")
  ITN.r <- rast(filename, win = ext(wgs.r))
  pITN <- project(ITN.r, crs(active.r))
  ITN.r <- resample(pITN, active.r, method = 'average')
  writeRaster(ITN.r, filename = paste0("ITN_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff")
}

