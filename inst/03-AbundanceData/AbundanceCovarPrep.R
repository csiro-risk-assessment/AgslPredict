# abundance aggregated to grid cells
# generanted from AbundancePrep.R
abund <- readRDS("abundance.rds")

# after generating abundance.rds
# get meteorological data downloaded from NASA POWER
# in directory meteo_for_abundance
# and spatial covariates
# in directory covariates_spatial
# to compile dataset for analysis of abundance

# weekly running means of precip -----------------------------------------------
# lagged from one day to duration of immature life history stage
# for each observation

# lat/long for each cell
precip.cells <- unique(abund$cell)
precip.lon <- precip.lat <- numeric(length(precip.cells))
for (i in 1:length(precip.cells)) {
  precip.lon[i] <- unique(abund[abund$cell == precip.cells[i], "lon.centroid"])
  precip.lat[i] <- unique(abund[abund$cell == precip.cells[i], "lat.centroid"])
}
precip.locations <- as.data.frame(cbind(precip.cells, precip.lon, precip.lat))

imm.duration <- 10 # duration of immature phase
p.duration <- 7 # duration of precip running mean

# running weekly averages
# for lag of one day to lag of imm.duration days wrt every observation
lag.p <- lag.t <- lag.r <- as.data.frame(matrix(NA, nrow = nrow(abund), ncol = imm.duration,
                                                dimnames = list(NULL, paste0("lag", 1:imm.duration))
))
names(lag.p) <- paste0("p_", names(lag.p))
names(lag.t) <- paste0("t_", names(lag.t))
names(lag.r) <- paste0("r_", names(lag.r))

pabund <- cbind(abund, lag.p, lag.t, lag.r)

## precipitations
for (i in 1:nrow(precip.locations)) {

  rec <- read.csv(paste0("../meteo_for_abundance/precip_", precip.cells[i], ".csv"),
                  header = TRUE)

  # subset dataframe to correspond with cell
  sub.cell <- pabund[pabund$cell == precip.cells[i], ]

  # lagged weekly averages of precip
  for (j in 1:nrow(sub.cell)) {
    id <- which(rec$YYYYMMDD == sub.cell$date[j])
    # lags from one day to imm.duration days
    for (k in 1:imm.duration) {
      # p.duration running means over lags of one day to imm.duration days
      sub.cell[j, colnames(lag.p)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "PRECTOTCORR"])
      sub.cell[j, colnames(lag.t)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "T2M"])
      sub.cell[j, colnames(lag.r)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "RH2M"])
    }
  }

  pabund[pabund$cell == precip.cells[i], ] <- sub.cell

}

all.equal(abund, pabund[ , names(abund)])

# spatial covariates -----------------------------------------------------------

# following uses same code as "relativeAbundancePrep.R"
# switching data.frame ra for abund

# Lambert azimuthal equal area
# see Steinwand et al. (1995)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

# use  limits denoted in first line of csv files
# should read: #xmin=-4099134.0	ymin=-4202349.0	cell_size=5000.0	nx=1520	ny=1280

# check active cells for correct extent info
all(read.csv("../covariates_spatial/active.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# read in spatial scope of RA in projected CRS
active <- read.csv("../covariates_spatial/active.csv", skip = 1, header = FALSE)
all.equal(dim(active), c(1280, 1520)) # expecting same number of rows and columns for all covariate data
active <- as.matrix(active)
library(terra)
# extent and crs
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

# read in elevation
elev.r <- rast("../covariates_spatial/elev.tif")

# read in covariates with existing .csv files (10 August 2023)

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

# assign covariate values to abundance data ------------------------------------

# project coords from abundance data and find corresponding cells
# (could also use cell id but sticking with geographic overlay)
coords <- cbind(abund$lon.centroid, abund$lat.centroid)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, crs)

plot(active.r)
plot(proj.coords, add = TRUE)

o.d2c <- extract(d2c.r, proj.coords)
o.d2r <- extract(d2r.r, proj.coords)
o.pop <- extract(pop.r, proj.coords) # no missing values
o.elev <- extract(elev.r, proj.coords) # no missing values

abund$d2c <- o.d2c[ , 2]
abund$d2r <- o.d2r[ , 2]
abund$pop <- o.pop[ , 2]
abund$elev <- o.elev[ , 2]

# bind spatial and temporal covariates -----------------------------------------

abundance <- cbind(abund, pabund[ , c(names(lag.p), names(lag.t), names(lag.r))])
saveRDS(abundance, file = "abund_obs_covs.rds")
