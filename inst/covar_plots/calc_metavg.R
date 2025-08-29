# spatial covariates
# Lambert azimuthal equal area
# see Steinwand et al. (1995)

library(terra)
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
# read active grid
active.r <- rast("../covariates_spatial/activeAfrica.tif")
res(active.r) # 5000 5000 with skip = 1

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
results.met <- cbind(
  mapCells,
  matrix(NA, nrow = nrow(mapCells), ncol = 29,
         dimnames = list(
           NULL,
           c(
             "lon", "lat",
             paste("T2M_", "Q", 1:4, sep = ""),
             paste("RH2M_", "Q", 1:4, sep = ""),
             paste("PRECTOTCOR_", "Q", 1:4, sep = ""),
             paste0(rep(c("T2M", "RH2M", "PRECTOTCOR"), each = 5),
                   c("_mean", paste0("_y", seq(2002, 2020, 6))))
           )
         ))
)

iter <- 1
N <- nrow(ap.crds)
pt <- proc.time()["elapsed"]

# for each MERRA-2 cell
for (u in u.cells) {

  ## meteorological covariates -------------------------------------------------
  met <- readRDS(paste0("../meteo_prediction_raw/met_",
                        m2crds[m2crds[ , "cell"] == u, "download.id"],
                        ".rds"))
  
  # prediction period is 2002--2020
  met <- met[215:7154, ]

  # identify quarters and years
  dates <- as.character(met[ , "YYYYMMDD"])
  quart <- lubridate::quarter(dates)
  years <- lubridate::year(dates)

  # for each active cell in a MERRA-2 cell
  au.cells <- results.met[results.met[ , "M2CellNum"] == u, "cell"]

  for (a in au.cells) {

    ## spatial covariates ------------------------------------------------------
    scovs <- c(
      pa.crds[a, ] # lon, lat
    )
    names(scovs) <- c("lon.centroid", "lat.centroid")

    # warning check on lon/lat
    if (abs(scovs["lon.centroid"] - met[1, "LON"]) > 5/8)
      stop("longitude for active cell unexpectedly far from MERRA-2 gridpoint")
    if (abs(scovs["lat.centroid"] - met[1, "LAT"]) > 0.5)
      stop("latitude for active cell unexpectedly far from MERRA-2 gridpoint")

      # assign results.met
      results.met[which(results.met[ , "cell"] == a), 5:33] <- c(
        # lon lat
        pa.crds[a, ],
        # quarterly means
        mean(met[quart == 1, "T2M"]),
        mean(met[quart == 2, "T2M"]),
        mean(met[quart == 3, "T2M"]),
        mean(met[quart == 4, "T2M"]),
        mean(met[quart == 1, "RH2M"]),
        mean(met[quart == 2, "RH2M"]),
        mean(met[quart == 3, "RH2M"]),
        mean(met[quart == 4, "RH2M"]),
        mean(met[quart == 1, "PRECTOTCORR"]),
        mean(met[quart == 2, "PRECTOTCORR"]),
        mean(met[quart == 3, "PRECTOTCORR"]),
        mean(met[quart == 4, "PRECTOTCORR"]),
        # overall means and annual means
        mean(met[ , "T2M"]),
        mean(met[years == 2002, "T2M"]),
        mean(met[years == 2008, "T2M"]),
        mean(met[years == 2014, "T2M"]),
        mean(met[years == 2020, "T2M"]),
        mean(met[ , "RH2M"]),
        mean(met[years == 2002, "RH2M"]),
        mean(met[years == 2008, "RH2M"]),
        mean(met[years == 2014, "RH2M"]),
        mean(met[years == 2020, "RH2M"]),
        mean(met[ , "PRECTOTCORR"]),
        mean(met[years == 2002, "PRECTOTCORR"]),
        mean(met[years == 2008, "PRECTOTCORR"]),
        mean(met[years == 2014, "PRECTOTCORR"]),
        mean(met[years == 2020, "PRECTOTCORR"])
      )
    }

    iter <- iter + 1
}

saveRDS(results.met, "results.met.rds")


