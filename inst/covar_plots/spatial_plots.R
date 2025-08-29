# This script reads the TIFF files of the following covariates:
# Elevation, distance to coast, probability of freshwater
# This script plots the map corresponding to each spatial covariate and saves
# the plots in the files elev.png, d2c.png and lakesrivers.png.

library(terra)
library(rnaturalearth)

# shoreline from rnaturalearth
shore <- ne_coastline(scale = "medium", returnclass = "sv")
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# read active grid
active.r <- rast("../covariates_spatial/activeAfrica.tif")
res(active.r)

# read in spatial covariate files
# elevation
elev.r <- rast("../covariates_spatial/elev.tif")
# distance to coast
d2c.r <- rast("../covariates_spatial/d2c.tif")
# probability of freshwater
lakesrivers.r <- rast("../covariates_spatial/lakesrivers.tif")

# check the resolution
res(d2c.r); res(lakesrivers.r); res(elev.r)

# distance to coast in km
d2c.r <- d2c.r/1000
# remove inactive cells
d2c.r[active.r==0] <- NA
lakesrivers.r[active.r==0] <- NA
elev.r[active.r==0] <- NA

# projection to WGS84
d2c.r <- project(d2c.r, "epsg:4326")
lakesrivers.r <- project(lakesrivers.r, "epsg:4326")
elev.r <- project(elev.r, "epsg:4326")

# plots ------------------------------------------------------------------------

library(rasterVis)
library(marmap)
library(latticeExtra)

brks.d2c <- c(200*(0:10))
tcks.d2c <- 0:10
labs.d2c <- c(200*(0:10))

brks.lakesrivers <- 0.1*c(0:10)
tcks.lakesrivers <- 0.1*c(0:10)
labs.lakesrivers <- 0.1*c(0:10)


brks.elev <- 500*c(0:8)
tcks.elev <- 500*c(0:8)
labs.elev <- c(expression("0"),expression("500"),expression("1000"),
               expression("1500"),expression("2500"),expression("2000"),
               expression("3000"),expression("3500"),expression("4000"))

cols <- c(RColorBrewer::brewer.pal(9, "GnBu")[9:1])

ann.theme.d2c <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)


png("d2c.png", 5*240, 5*240*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(d2c.r,
          main =  "Distance to nearest coast",
          par.settings = ann.theme.d2c,
          at = brks.d2c,
          colorkey = list(
            at = tcks.d2c,
            labels = list(
              at = tcks.d2c,
              labels = labs.d2c
            )
          ),
          margin = FALSE
) + layer(llines(shore, col = "black"))
dev.off()


cols <- c(RColorBrewer::brewer.pal(9, "BrBG"))

ann.theme.prob <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

png("lakesrivers.png", 5*240, 5*240*(0.8),pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(lakesrivers.r,
          main =  "Probability of permanent freshwater",
          par.settings = ann.theme.prob,
          at = brks.lakesrivers,
          colorkey = list(
            at = tcks.lakesrivers,
            labels = list(
              at = tcks.lakesrivers,
              labels = labs.lakesrivers
            )
          ),
          margin = FALSE
) + layer(llines(shore, col = "black"))
dev.off()


cols <- rev(etopo.colors(11))[c(3:10)]
ann.theme.elev <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)
png("elev.png", 5*240, 5*240*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(elev.r,
          main =  "Bedrock elevation",
          par.settings = ann.theme.elev,
          at = brks.elev,
          colorkey = list(
            at = tcks.elev,
            labels = list(
              at = tcks.elev,
              labels = labs.elev
            )
          ),
          margin = FALSE
) + layer(llines(shore, col = "black"))
dev.off()
