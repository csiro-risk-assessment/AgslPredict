# outputs:
# for each 5 km x 5 km cell
# predictions per (financial) quarter for:
#   - species abundances
# predictions averaged per year:
#   - Pianka proportional resource usage terms
# for a subset of cells:
#   - full time series of daily abundances

# inputs:
#   - multinomial fit
#   - spatial covariates for multinomial
#   - met covariates for multinomial
#   - abundance model parameters
#   - spatial covariates for abundance
#   - scalings for covariates
#   - precip for abundance
#   - prediction function pxPrediction

# for array job set id <- NULL
id <- NULL
if (is.null(id)) {
  args <- commandArgs(trailingOnly=T)
  id <- as.numeric(args[1])
}

# Read in multinomial fit, abundance parameters and spatial covariates ---------

# prediction function
source("pxPrediction_function.R")

# multinomial model fit
m <- readRDS("../02-RelativeAbundanceAnalysis/multinomialGLMfit1m.rds")
ra.scaling <- readRDS("../02-RelativeAbundanceAnalysis/raScaling.rds")
ra.pars <- readRDS("../02-RelativeAbundanceAnalysis/raPars.rds")

# abundance model fit and covariate scalings
a.mod <- readRDS("../04-AbundanceAnalysis/abundanceGLMfit.rds")
a.scaling <- readRDS("../04-AbundanceAnalysis/abundanceScaling.rds")

# abundance predictions correspond to indoor hlc
# observation covariates set to zero
obs.names <- c("light", "exit", "PSC", "animal",  "pit",  "outdoor")

# spatial covariates
# Lambert azimuthal equal area
# see Steinwand et al. (1995)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

library(terra)

#  vector intervention covariates (spatio-temporal): ITNs ----------------------

YITN <- 2002:2020
ITN.files <- paste0("../vector_intervention/ITN_", YITN, ".tif")
itn.r <- rast(ITN.files)
names(itn.r) <- paste0("Y", YITN)

active.r <- rast("../covariates_spatial/activeAfrica.tif")
same.crs(crs(active.r), crs) # TRUE
res(active.r)

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

# meteorological covariates ----------------------------------------------------

# predictions will use CRS from meteorological covariates (WGS84)

# meteorological covariate grid and download ids
# raster
m2.r <- readRDS("../meteo_prediction_raw/MERRA2raster.rds")

# raster of results ------------------------------------------------------------

# MERRA2 cell ids
a.crds <- crds(active.r)
a.crds <- cbind(a.crds, 1:nrow(a.crds))
colnames(a.crds)[3] <- "cell"
# subset to active cells
va <- ifelse(values(active.r) == 0, NA, values(active.r))
a.crds <- a.crds[!is.na(va), ]
# project to m2
ap.crds <- vect(as.data.frame(a.crds), geom = c("x", "y"), crs = crs(active.r))
ap.crds <- project(ap.crds, crs(m2.r))
# get cell number of MERRA-2 for each active centroid
am2cells <- extract(m2.r, ap.crds, cells = TRUE, ID = FALSE)
# mapping of cell numbers between active / spatial covars and MERRA-2
mapCells <- cbind(a.crds, am2cells$cell)
colnames(mapCells)[4] <- "M2CellNum"

# load download ids for raw meteorological from NASA POWER
m2crds <- read.csv("../meteo_prediction_raw/MERRA-2-coordinates.csv", header = TRUE)

# prediction -------------------------------------------------------------------

# for each active cell in WGS84
# find subset with shared meteorological covariates
# read in met covars and process
# get spatial covariates for each active cell with shared met covars
# calculate relative abundance, abundance and Pianka components
# export using active cells in WGS84

# time period for prediction is 1 Jan 2002 to 31 Dec 2020

# MERRA-2 cells that include active cells
u.cells <- sort(unique(mapCells[ , "M2CellNum"]))

# for checking:
check.crds <- crds(active.r)
pact.v <- vect(check.crds, crs = crs(active.r))
proj.a <- project(pact.v, crs(m2.r))
pa.crds <- crds(proj.a)

# random sample for full time series
set.seed(1000)
full.ids <- sample(mapCells[ , "cell"], 500)

# store results
results <- cbind(
  mapCells,
  matrix(NA, nrow = nrow(mapCells), ncol = 36,
         dimnames = list(
           NULL,
           c(
             "lon", "lat",
             paste("Aa_", "Q", 1:4, sep = ""),
             paste("Ac_", "Q", 1:4, sep = ""),
             paste("Ag_", "Q", 1:4, sep = ""),
             paste0(rep(c("Aa", "Ac", "Ag"), each = 5),
                   c("_mean", paste0("_y", seq(2002, 2020, 6)))),
             "Aa_r", "Ac_r", "Ag_r",
             "logAcAa_v", "logAgAa_v",
             "tot_mean", "tot_v"
           )
         ))
)

iter <- 1
N <- nrow(ap.crds)
pt <- proc.time()["elapsed"]

start.id <- ((id - 1)*1000 + 1)
end.id <- min(id*1000, length(u.cells))

u.cells <- u.cells[start.id:end.id]

# for each MERRA-2 cell
for (u in u.cells) {

  ## meteorological covariates -------------------------------------------------
  met <- readRDS(paste0("../meteo_prediction_raw/met_",
                        m2crds[m2crds[ , "cell"] == u, "download.id"],
                        ".rds"))

  ### abundance uses 10 day (immature duration) average of
  # weekly running means of precip
  # row 215 of met is 1 Jan 2002
  # row 7154 of met is 31 Dec 2020
  immprecip <- cbind(
    filter(met[(215 - 7):(7154 - 1), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 1
    filter(met[(215 - 8):(7154 - 2), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 2
    filter(met[(215 - 9):(7154 - 3), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 3
    filter(met[(215 - 10):(7154 - 4), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 4
    filter(met[(215 - 11):(7154 - 5), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 5
    filter(met[(215 - 12):(7154 - 6), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 6
    filter(met[(215 - 13):(7154 - 7), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 7
    filter(met[(215 - 14):(7154 - 8), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 8
    filter(met[(215 - 15):(7154 - 9), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1), # lag 9
    filter(met[(215 - 16):(7154 - 10), c("T2M", "PRECTOTCORR")],
           rep(1/7, 7), sides = 1) # lag 10
  )
  immprecip <- immprecip[7:nrow(immprecip), ] # drop first 7 days of running means
  colnames(immprecip) <- paste0(c("t_lag", "p_lag"), rep(1:10, each = 2))

  pmeans <- rowMeans(immprecip[ , grep("p_lag", colnames(immprecip))])
  tmeans <- rowMeans(immprecip[ , grep("t_lag", colnames(immprecip))])

  # observation covariates: known zero values
  obs.cov <- matrix(0, nrow = nrow(immprecip), ncol = length(obs.names),
                    dimnames = list(NULL, obs.names))

  ### relative abundance 30 day running averages
  ra.1m <- stats::filter(met[186:7154, c("T2M", "RH2M", "PRECTOTCORR")],
                         rep(1/30, 30), sides = 1)
  rownames(ra.1m) <- as.character(met[186:7154, "YYYYMMDD"])
  colnames(ra.1m) <- paste0(c("T2M", "RH2M", "PRECTOTCORR"), "1m")
  ra.1m <- ra.1m[30:nrow(ra.1m), ]

  # identify quarters and years
  dates <- rownames(ra.1m)
  quart <- lubridate::quarter(dates)
  years <- lubridate::year(dates)

  # for each active cell in a MERRA-2 cell
  au.cells <- results[results[ , "M2CellNum"] == u, "cell"]

  for (a in au.cells) {

    ## spatial covariates ------------------------------------------------------
    scovs <- c(
      as.numeric(d2c.r[a]),
      as.numeric(pop.r[a]),
      as.numeric(d2r.r[a]),
      as.numeric(elev.r[a]),
      pa.crds[a, ] # lon, lat
    )
    names(scovs) <- c("d2c", "pop", "d2r", "elev", "lon.centroid", "lat.centroid")

    # warning check on lon/lat
    if (abs(scovs["lon.centroid"] - met[1, "LON"]) > 5/8)
      stop("longitude for active cell unexpectedly far from MERRA-2 gridpoint")
    if (abs(scovs["lat.centroid"] - met[1, "LAT"]) > 0.5)
      stop("latitude for active cell unexpectedly far from MERRA-2 gridpoint")

    if (!any(is.na(scovs)) && !is.na(as.numeric(itn.r[[1]][a]))) {

      scovs["lat.centroid"] <- abs(scovs["lat.centroid"])
      scovs.mat <- matrix(scovs, byrow = TRUE, nrow = nrow(ra.1m), ncol = length(scovs),
                          dimnames = list(
                            rownames(ra.1m), names(scovs)
                          ))

      # ITN data by year
      itn.df <- data.frame(YITN = YITN, ITN = as.numeric(itn.r[a]))
      ITN <- rep(NA, nrow(scovs.mat))
      scovs.mat <- cbind(scovs.mat, ITN)
      for (k in 1:length(years)) {
        scovs.mat[years == YITN[k], "ITN"] <- itn.df[k, "ITN"]
      }

      # relative abundance prediction data
      ra.dat <- as.data.frame(cbind(
        scovs.mat,
        ra.1m
      ))
      # origin year = 2020
      ra.dat$origin.year <- 2020
      # scale
      ra.names <- colnames(ra.dat)
      for (i in ra.names) {
        ra.dat[ , i] <- (ra.dat[ , i] - ra.scaling$mins[i])/(ra.scaling$maxs[i] - ra.scaling$mins[i])
      }
      ra.df <- as.data.frame(ra.dat)

      # abundance spatial prediction data
      a.mat <- scovs.mat[ , c("d2c", "d2r", "pop", "elev", "lat.centroid", "lon.centroid", "ITN")]
      # relabel lat/long, itn
      colnames(a.mat)[which(colnames(a.mat) == "lat.centroid")] <- "lat"
      colnames(a.mat)[which(colnames(a.mat) == "lon.centroid")] <- "long"
      colnames(a.mat)[which(colnames(a.mat) == "ITN")] <- "itn"
      # abs value of lat
      a.mat[ , "lat"] <- abs(a.mat[ , "lat"])
      # precip for abundance
      a.mat <- cbind(a.mat, pmeans)
      colnames(a.mat)[ncol(a.mat)] <- "precip"
      a.mat <- cbind(a.mat, tmeans)
      colnames(a.mat)[ncol(a.mat)] <- "temp"
      # scale
      a.names <- colnames(a.mat)
      for (i in a.names) {
        a.mat[ , i] <- (a.mat[ , i] - a.scaling$mins[i])/(a.scaling$maxs[i] - a.scaling$mins[i])
      }
      # bind observations covariates
      # offset "units" is equal to one
      a.df <- as.data.frame(cbind(obs.cov, a.mat))
      a.df$units <- 1

      # predictions:
      preds.ls <- pxPrediction(ra.df, a.df, m, ra.pars, a.mod, scovs["pop"])

      # actual abundance
      x <- preds.ls$x

      # assign results
      results[which(results[ , "cell"] == a), 5:40] <- c(
        # lon lat
        pa.crds[a, ],
        # quarterly means
        mean(x[quart == 1, "Aa"]),
        mean(x[quart == 2, "Aa"]),
        mean(x[quart == 3, "Aa"]),
        mean(x[quart == 4, "Aa"]),
        mean(x[quart == 1, "Ac"]),
        mean(x[quart == 2, "Ac"]),
        mean(x[quart == 3, "Ac"]),
        mean(x[quart == 4, "Ac"]),
        mean(x[quart == 1, "Ag"]),
        mean(x[quart == 2, "Ag"]),
        mean(x[quart == 3, "Ag"]),
        mean(x[quart == 4, "Ag"]),
        # overall means and annual means
        mean(x[ , "Aa"]),
        mean(x[years == 2002, "Aa"]),
        mean(x[years == 2008, "Aa"]),
        mean(x[years == 2014, "Aa"]),
        mean(x[years == 2020, "Aa"]),
        mean(x[ , "Ac"]),
        mean(x[years == 2002, "Ac"]),
        mean(x[years == 2008, "Ac"]),
        mean(x[years == 2014, "Ac"]),
        mean(x[years == 2020, "Ac"]),
        mean(x[ , "Ag"]),
        mean(x[years == 2002, "Ag"]),
        mean(x[years == 2008, "Ag"]),
        mean(x[years == 2014, "Ag"]),
        mean(x[years == 2020, "Ag"]),
        # ra raw predictions
        preds.ls$r,
        # ra variances
        preds.ls$r.v,
        # tot mean
        preds.ls$xi.m,
        # tot variance
        preds.ls$xi.v
      )

      if (a%in%full.ids) {
        saveRDS(x, file = paste0("x_", a, ".rds"))
      }

    }

    if (iter%%10 == 0)
      cat("iteration ", iter, " of ", N, ". Elapsed time: ",
          proc.time()["elapsed"] - pt, "\n", sep = "")

    iter <- iter + 1

  }
}

saveRDS(results, paste0("results_", id, ".rds"))

