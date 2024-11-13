# spatial covariates -----------------------------------------------------------

# spatial covariates
# Lambert azimuthal equal area
# see Steinwand et al. (1995)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# check active cells for correct extent info
all(read.csv("../covariates_spatial/active.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# read in spatial scope of RA in projected CRS
active <- read.csv("../covariates_spatial/active.csv", skip = 1, header = FALSE)
all.equal(dim(active), c(1280, 1520)) # expecting same number of rows and columns for all covariate data
active <- as.matrix(active)
library(terra)
# extent and crs from Nick's R file above
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

## read in elevation

elev.r <- rast("../covariates_spatial/elev.tif")

## read in remaining spatial covariate .csv files

# check for correct extent info
all(read.csv("../covariates_spatial/d2c.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
all(read.csv("../covariates_spatial/d2r.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
all(read.csv("../covariates_spatial/pop.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))

d2c <- read.csv("../covariates_spatial/d2c.csv", skip = 1, header = FALSE)
d2r <- read.csv("../covariates_spatial/d2r.csv", skip = 1, header = FALSE)
pop <- read.csv("../covariates_spatial/pop.csv", skip = 1, header = FALSE)

all.equal(dim(d2c), c(1280, 1520))
all.equal(dim(d2r), c(1280, 1520))
all.equal(dim(pop), c(1280, 1520))

d2c <- as.matrix(d2c)
d2r <- as.matrix(d2r)
pop <- as.matrix(pop)

d2c.r <- rast(d2c[nrow(d2c):1, ], crs = crs,
              extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
d2r.r <- rast(d2r[nrow(d2r):1, ], crs = crs,
              extent =  c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
pop.r <- rast(pop[nrow(pop):1, ], crs = crs,
              extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))

res(d2c.r); res(d2r.r); res(pop.r); res(elev.r)


d2c.r[active.r==0] <- NA
d2r.r[active.r==0] <- NA
pop.r[active.r==0] <- NA
elev.r[active.r==0] <- NA


# projection to WGS84
d2c.r <- project(d2c.r, "epsg:4326")
d2r.r <- project(d2r.r, "epsg:4326")
pop.r <- project(pop.r, "epsg:4326")
elev.r <- project(elev.r, "epsg:4326")

# plots ------------------------------------------------------------------------

library(rasterVis)

RColorBrewer::display.brewer.all(n=7, type="seq", select=NULL, exact.n=TRUE,
                                                     colorblindFriendly=TRUE)
RColorBrewer::display.brewer.all(n=7, type="div", select=NULL, exact.n=TRUE,
                                 colorblindFriendly=TRUE)

## annual ----------------------------------------------------------------------

brks.d2c <- c(200*(0:10))
tcks.d2c <- 0:10
labs.d2c <- c(200*(0:10))

brks.d2r <- 0.5*c(1:12)
tcks.d2r <- 0.5*c(1:12)
labs.d2r <- c("", expression(10^1), "", expression(10^2),"",
              expression(10^3), "", expression(10^4), "",
              expression(10^5), "", expression(10^6))


brks.elev <- 500*c(0:8)
tcks.elev <- 500*c(0:8)
labs.elev <- c(expression("0"),expression("500"),expression("1000"),
               expression("1500"),expression("2500"),expression("2000"),
               expression("3000"),expression("3500"),expression("4000"))


brks.pop <- 0.5*c(0:10)
tcks.pop <- 0.5*c(0:10)
labs.pop <- c("", 5, 10, 50, 100, 500, 1000, 5000,
              10000, 50000, ""
)

cols <- c(RColorBrewer::brewer.pal(9, "GnBu")[9:1])

ann.theme.d2c <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

cols <- c(RColorBrewer::brewer.pal(11, "RdYlGn")[11:1])
ann.theme.pop <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

library(marmap)

cols <- rev(etopo.colors(11))[c(3:10)]
ann.theme.elev <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

png("d2c.png", 2.5*240, 2.5*240, pointsize = 72)
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
)
dev.off()


png("d2r.png", 2.5*240, 2.5*240, pointsize = 72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(log10(d2r.r),
          main =  "Distance to nearest river",
          par.settings = ann.theme.d2c,
          at = brks.d2r,
          colorkey = list(
            at = tcks.d2r,
            labels = list(
              at = tcks.d2r,
              labels = labs.d2r
            )
          ),
          margin = FALSE
)
dev.off()


pop.rlog <- log10(pop.r)
pop.rlog[pop.r<=1] <- 0
png("pop.png", 2.5*240, 2.5*240, pointsize = 72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(pop.rlog,
          main =  "Human population density",
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
)
dev.off()


png("elev.png", 2.5*240, 2.5*240, pointsize = 72)
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
)
dev.off()
