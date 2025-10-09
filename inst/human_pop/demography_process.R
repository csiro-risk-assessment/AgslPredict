library(terra)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
active.r <- rast("../covariates_spatial/activeAfrica.tif")

# project to WGS84
wgs.r <- project(active.r, "epsg:4326")
print(ext(active.r))
active.mask <- active.r != 0


# download raw files into preprocessing directory
# available from zip files
#
# "landscan-global-2000.tif" <https://doi.org/10.48690/1524196>
# "landscan-global-2001.tif" <https://doi.org/10.48690/1524197>
# "landscan-global-2002.tif" <https://doi.org/10.48690/1524198>
# "landscan-global-2003.tif" <https://doi.org/10.48690/1524199>
# "landscan-global-2004.tif" <https://doi.org/10.48690/1524200>
# "landscan-global-2005.tif" <https://doi.org/10.48690/1524201>
# "landscan-global-2006.tif" <https://doi.org/10.48690/1524202>
# "landscan-global-2007.tif" <https://doi.org/10.48690/1524203>
# "landscan-global-2008.tif" <https://doi.org/10.48690/1524204>
# "landscan-global-2009.tif" <https://doi.org/10.48690/1524205>
# "landscan-global-2010.tif" <https://doi.org/10.48690/1524206>
# "landscan-global-2011.tif" <https://doi.org/10.48690/1524207>
# "landscan-global-2012.tif" <https://doi.org/10.48690/1524215>
# "landscan-global-2013.tif" <https://doi.org/10.48690/1524208>
# "landscan-global-2014.tif" <https://doi.org/10.48690/1524209>
# "landscan-global-2015.tif" <https://doi.org/10.48690/1524210>
# "landscan-global-2016.tif" <https://doi.org/10.48690/1524211>
# "landscan-global-2017.tif" <https://doi.org/10.48690/1524212>
# "landscan-global-2018.tif" <https://doi.org/10.48690/1524213>
# "landscan-global-2019.tif" <https://doi.org/10.48690/1524214>
# "landscan-global-2020.tif" <https://doi.org/10.48690/1523378>
# "landscan-global-2021.tif" <https://doi.org/10.48690/1527702>

preprocessing.dir <- getwd()  # for example if downloaded and unzipped into working directory

#### To be used on High Peformance Computing (HPC): use the number associated with process ####

# Example with Slurm:
# SLURM_ARRAY_TASK_ID = Sys.getenv("SLURM_ARRAY_TASK_ID") # for HPC parallel control
# SLURM_ARRAY_TASK_ID = as.numeric(SLURM_ARRAY_TASK_ID) # making sure slurm_array_id is numeric
# # argument number between 1 and 22 (the number of  years from 2000 to 2021)
# no.batch = SLURM_ARRAY_TASK_ID

##### If run on laptop or PC process each file in a separate session
##### where no.batch varies from 1 to 22
no.batch <- 1 # no.batch varies from 1 to 22

# getting the year and the file to process
year <- no.batch - 1 + 2000
f <- paste0(preprocessing.dir, "/landscan-global-", year, ".tif")
print(paste0("The year of the human population to process is ", year, " from the file : ", f))

#read the file in a rast
dem.r <- rast(f, win = ext(wgs.r))

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

####  clean up
rm(dem.r, dem.proj, dem.resamp.sum, dem.proj.size, dem.size)
gc()
####

# write new file
writeRaster(dem.sum.r, filename = paste0("./human_pop_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff", overwrite = TRUE)

#### clean up
rm(dem.sum.r)
gc()
Sys.sleep(15)
####

# write new file
writeRaster(dem.avg.r, filename = paste0("./human_density_", year, ".tif"), datatype = "FLT8S", filetype = "GTiff", overwrite = TRUE)
