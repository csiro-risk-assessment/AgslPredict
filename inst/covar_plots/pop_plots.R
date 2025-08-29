# This script reads the TIFF files of the human population density for years
# 2002, 2008, 2014 and 2020, plots the maps corresponding to human density on
# a 2 by 2 layout and saves the plots in the file Yearly_pop_density.png

library(terra)
library(rnaturalearth)

# shoreline from rnaturalearth
shore <- ne_coastline(scale = "medium", returnclass = "sv")
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# read active grid
active.r <- rast("../covariates_spatial/activeAfrica.tif")
res(active.r)

# raster layers
pop.r <- rast(active.r,
                  nlyrs = 4,
                  names = c(
                    paste0(rep(c("human_density"), each = 4),
                           c(paste0("_", seq(2002, 2020, 6))))
                  )
)

poptif.dir <- "../human_pop/"
years = seq(2002, 2020, 6)

for (year in years)
{
  pop.r[[paste0("human_density_",year)]] <- rast(paste0("../human_pop/human_density_", year, ".tif"))
}

pop.r <- project(pop.r, "epsg:4326")

# plots ------------------------------------------------------------------------

library(rasterVis)
library(marmap)
library(latticeExtra)

brks.pop <- 0.5*c(0:10)
tcks.pop <- 0.5*c(0:10)
labs.pop <- c("", 5, 10, 50, 100, 500, 1000, 5000,
              10000, 50000, ""
)
cols <- c(RColorBrewer::brewer.pal(11, "RdYlGn")[11:1])
ann.theme.pop <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)

pop.rlog <- log10(pop.r)
pop.rlog[pop.r<=1] <- 0

png("Yearly_pop_density.png", 5*240*3, 5*240*3*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(pop.rlog,
          main =  "Human population density",
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = ann.theme.pop,
          at = brks.pop,
          colorkey = list(
            at = tcks.pop,
            labels = list(
              at = tcks.pop,
              labels = labs.pop
            )
          ),
          margin = FALSE
) + layer(llines(shore, col = "black"))
dev.off()
