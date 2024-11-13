# outputs:
# for each meteorogloical cell
# predictions of avg length of durations with temperatures
#   - 18C or below
#   - 34C or above

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
# extent and crs
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

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
# export using active cells in WGS84

# time period for prediction is 1 Jan 2002 to 31 Dec 2020

# MERRA-2 cells that include active cells
u.cells <- sort(unique(mapCells[ , "M2CellNum"]))

# for checking:
check.crds <- crds(active.r)
pact.v <- vect(check.crds, crs = crs(active.r))
proj.a <- project(pact.v, crs(m2.r))
pa.crds <- crds(proj.a)

# store results
resultsT <- cbind(
  mapCells,
  matrix(NA, nrow = nrow(mapCells), ncol = 1,
         dimnames = list(
           NULL,
           c(
             "TOK"
           )
         ))
)

iter <- 1
N <- nrow(ap.crds)
pt <- proc.time()["elapsed"]

m2.r$Ttol <- NA

# for each MERRA-2 cell
for (u in u.cells) {

  ## meteorological covariates -------------------------------------------------
  met <- readRDS(paste0("../meteo_prediction_raw/met_",
                        m2crds[m2crds[ , "cell"] == u, "download.id"],
                        ".rds"))

  ### for 1 Jan 2002 to 31 Dec 2020
  # calculate avg interval length where daily either
  # is either 34C or above
  # or 18C or below
  metsub <- met[met$YEAR > 2001, ]
  metsub <- metsub[metsub$YEAR < 2021, ]

  T.OK <- metsub$T2M >= 18 & metsub$T2M <= 34
  dNOK <- filter(!T.OK, rep(1, 30), sides = 1)

  # # if half years have interval of temp tol exceedance
  # # then no pass
  # avgT <- mean(tapply((dNOK == 30), list(metsub$YEAR), sum, na.rm = TRUE) > 0)
  # m2.r$Ttol[u] <- avgT >= 0.50
  # if every year has at least one 30 day interval of temp tol exceedance
  # then no pass, pass threshold otherwise
  pass <- !all(tapply((dNOK == 30), list(metsub$YEAR), sum, na.rm = TRUE) > 0)
  m2.r$Ttol[u] <- pass

  # for each active cell in a MERRA-2 cell
  resultsT[resultsT[ , "M2CellNum"] == u, "TOK"] <- pass

  if (iter%%10 == 0)
    cat("iteration ", iter, " of ", length(u.cells), ". Elapsed time: ",
        proc.time()["elapsed"] - pt, "\n", sep = "")

  iter <- iter + 1

}

saveRDS(resultsT, "threshmat.rds")
saveRDS(m2.r, file = "Tthreshold.rds")


