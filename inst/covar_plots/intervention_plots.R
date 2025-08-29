# This script reads the TIFF files of ITN and IRS for years
# 2002, 2008, 2014 and 2020, plots the maps corresponding to ITN and IRS on
# two different 2 by 2 layouts and saves the plots in the files
# Yearly_ITN.png and Yearly_IRS.png

library(terra)
library(rnaturalearth)

# shoreline from rnaturalearth
shore <- ne_coastline(scale = "medium", returnclass = "sv")
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# read active grid
active.r <- rast("../covariates_spatial/activeAfrica.tif")
res(active.r)

# raster layers
intervention.r <- rast(active.r,
              nlyrs = 8,
              names = c(
                paste0(rep(c("ITN", "IRS"), each = 4),
                       c(paste0("_", seq(2002, 2020, 6))))
              )
)

intervention.dir <- "../vector_intervention/"
years = seq(2002, 2020, 6)

for (year in years)
{
  intervention.r[[paste0("ITN_",year)]] <- rast(paste0(intervention.dir, "ITN_", year, ".tif"))
  intervention.r[[paste0("IRS_",year)]] <- rast(paste0(intervention.dir, "IRS_", year, ".tif"))
}

# projection to WGS84
intervention.r <- project(intervention.r, "epsg:4326")

# plots ------------------------------------------------------------------------
library(rasterVis)
library(marmap)
library(latticeExtra)

cols <- rev(c(RColorBrewer::brewer.pal(9, "GnBu")[9:1]))
ann.theme.itn <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)
cols <- rev(c(RColorBrewer::brewer.pal(9, "BuPu")[9:1]))
ann.theme.irs <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

png("Yearly_ITN.png", 2.5*480*3, 2.5*480*3*(0.8), pointsize = 72, res = 6*72)
levelplot(intervention.r[[c("ITN_2002","ITN_2008","ITN_2014","ITN_2020")]],
          main = expression("ITN use mean"),
          names.attr = paste0("Year ", c("2002","2008","2014","2020")),
          layout = c(2, 2),
          par.settings = ann.theme.itn
)  + layer(llines(shore, col = "black"))
dev.off()




png("Yearly_IRS.png", 2.5*480*3, 2.5*480*3*(0.8), pointsize = 72, res = 6*72)
levelplot(intervention.r[[c("IRS_2002", "IRS_2008", "IRS_2014", "IRS_2020")]],
          main = expression("IRS use mean"), cex = 2,
          names.attr = paste0("Year ", c("2002","2008","2014","2020")),
          layout = c(2, 2),
          par.settings = ann.theme.irs
) + layer(llines(shore, col = "black"))
dev.off()

