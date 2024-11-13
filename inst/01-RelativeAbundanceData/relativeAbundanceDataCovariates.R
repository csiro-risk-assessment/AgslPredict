# VectorBase relative abundance data and covariates ----------------------------
# zipping together VectorBase relative abundance data
# compiled to grid cell resolution

# read in relative abundance data previously output to this directory
# by running "RelativeAbundanceDatePrep.R" also in this directory
ra <- readRDS("VB_PCR_rel_abundance_stra.rds")

# non meteorological covariates ------------------------------------------------

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

## read in elevation -----------------------------------------------------------

elev.r <- rast("../covariates_spatial/elev.tif")

## read in remaining spatial covariate .csv files ------------------------------

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

## assign covariate values to relative abundance data ---------------------------

# project coords from abundance data and find corresponding cells
# (could also use cell id but sticking with geographic overlay)
coords <- cbind(ra$lon.centroid, ra$lat.centroid)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, crs)
par(mfrow = c(1, 1))
plot(active.r)
plot(proj.coords, add = TRUE)

o.d2c <- extract(d2c.r, proj.coords)
sum(o.d2c[ , 2] >= 0) == nrow(proj.coords)
o.d2r <- extract(d2r.r, proj.coords)
sum(o.d2r[ , 2] >= 0) == nrow(proj.coords)
o.pop <- extract(pop.r, proj.coords)
sum(o.pop[ , 2] >= 0) == nrow(proj.coords)
o.elev <- extract(elev.r, proj.coords) # no missing values

ra$d2c <- o.d2c[ , 2]
ra$d2r <- o.d2r[ , 2]
ra$pop <- o.pop[ , 2]
ra$elev <- o.elev[ , 2]

# meteorological covariates ----------------------------------------------------

# test model: two datasets
# one with 31 day max duration for sample (~1 month)
# one with 93 day max duration for sample (~3 months)

ra3m <- ra
ra1m <- ra[ra$diffDate <= 31, ]
rm(ra)

ra.names <- c("ra3m", "ra1m")
ra.ls <- list(ra3m = NULL, ra1m = NULL)

for (r in ra.names) {

  ra <- get(r)

  # for each choice of max duration construct met covariates and append

  # lat/long for each cell
  meteo.cells <- unique(ra$cell)
  meteo.lon <- meteo.lat <- numeric(length(meteo.cells))
  for (i in 1:length(meteo.cells)) {
    meteo.lon[i] <- unique(ra[ra$cell == meteo.cells[i], "lon.centroid"])
    meteo.lat[i] <- unique(ra[ra$cell == meteo.cells[i], "lat.centroid"])
  }
  meteo.locations <- as.data.frame(cbind(meteo.cells, meteo.lon, meteo.lat))

  # we store the average of PRECTOTCORR and T2M and RH2M for the last month (30 days), last
  # three months (90 days) and the last year (360)
  avg.p <- as.data.frame(matrix(NA, nrow = nrow(ra), ncol = 6))
  colnames(avg.p) <- c("PRECTOTCORR1m", "T2M1m", "RH2M1m",
                       "PRECTOTCORR3m", "T2M3m", "RH2M3m")
  ra <- cbind(ra, avg.p)

  for (i in 1:nrow(meteo.locations)) {

    if (file.exists(paste0("../meteo_for_relative_abundance/meteo_", meteo.cells[i], ".csv"))) {
      rec <- read.csv(paste0("../meteo_for_relative_abundance/meteo_", meteo.cells[i], ".csv"),
                      header = TRUE)
      # subset dataframe to correspond with cell
      sub.cell <- ra[ra$cell == meteo.cells[i], ]
      # lagged monthly averages of precip
      for (j in 1:nrow(sub.cell)) {
        date1 <- sub.cell$date1[j]
        date2 <- sub.cell$date2[j]
        id1 <- which(rec$YYYYMMDD == as.Date(date1))
        id2 <- which(rec$YYYYMMDD == as.Date(date2))
        meanPr1m <- mean(rec[(id1 - 30 + 1):id1, "PRECTOTCORR"])
        meanT2M1m <- mean(rec[(id1 - 30 + 1):id1, "T2M"])
        meanRH2M1m <- mean(rec[(id1 - 30 + 1):id1, "RH2M"])
        meanPr3m <- mean(rec[(id1 - 90 + 1):id1, "PRECTOTCORR"])
        meanT2M3m <- mean(rec[(id1 - 90 + 1):id1, "T2M"])
        meanRH2M3m <- mean(rec[(id1 - 90 + 1):id1, "RH2M"])
        diffDate <- as.numeric(lubridate::date(date2) - lubridate::date(date1))
        if (diffDate > 1)
        {
          meanPr1m <- (mean(rec[(id2 - diffDate - 30 + 1):id2, "PRECTOTCORR"]) + meanPr1m)/2
          meanT2M1m <- (mean(rec[(id2 - diffDate - 30 + 1):id2, "T2M"]) + meanT2M1m)/2
          meanRH2M1m <- (mean(rec[(id2 - diffDate - 30 + 1):id2, "RH2M"]) + meanRH2M1m)/2
          meanPr3m <- (mean(rec[(id2 - diffDate - 90 + 1):id2, "PRECTOTCORR"]) + meanPr3m)/2
          meanT2M3m <- (mean(rec[(id2 - diffDate - 90 + 1):id2, "T2M"]) + meanT2M3m)/2
          meanRH2M3m <- (mean(rec[(id2 - diffDate - 90 + 1):id2, "RH2M"]) + meanRH2M3m)/2
        }
        sub.cell[j, "PRECTOTCORR1m"] <- meanPr1m
        sub.cell[j, "T2M1m"] <- meanT2M1m
        sub.cell[j, "RH2M1m"] <- meanRH2M1m
        sub.cell[j, "PRECTOTCORR3m"] <- meanPr3m
        sub.cell[j, "T2M3m"] <- meanT2M3m
        sub.cell[j, "RH2M3m"] <- meanRH2M3m
      }

      ra[ra$cell == meteo.cells[i], ] <- sub.cell

    } else {
      cat("missing: ", paste0("../meteo_for_relative_abundance/meteo_", meteo.cells[i], ".csv \n"))
    }

  }

  ra.ls[[r]] <- ra

}

par(mfrow = c(3, 2))
hist(ra.ls[["ra3m"]]$PRECTOTCORR1m)
hist(ra.ls[["ra3m"]]$T2M1m)
hist(ra.ls[["ra3m"]]$RH2M1m)
hist(ra.ls[["ra3m"]]$PRECTOTCORR3m)
hist(ra.ls[["ra3m"]]$T2M3m)
hist(ra.ls[["ra3m"]]$RH2M3m)

par(mfrow = c(3, 2))
hist(ra.ls[["ra1m"]]$PRECTOTCORR1m)
hist(ra.ls[["ra1m"]]$T2M1m)
hist(ra.ls[["ra1m"]]$RH2M1m)
hist(ra.ls[["ra1m"]]$PRECTOTCORR3m)
hist(ra.ls[["ra1m"]]$T2M3m)
hist(ra.ls[["ra1m"]]$RH2M3m)

# correlations between two choices of max duration
s <- intersect(ra.ls[["ra3m"]]$samp.id, ra.ls[["ra1m"]]$samp.id)
ra3mSub1m <- ra.ls[["ra3m"]][ra.ls[["ra3m"]]$samp.id%in%ra.ls[["ra1m"]]$samp.id, ]
names(ra3mSub1m) <- paste0(names(ra3mSub1m), "3")
pairs(cbind(
  ra3mSub1m[ , c("PRECTOTCORR1m3", "T2M1m3", "RH2M1m3",
                    "PRECTOTCORR3m3", "T2M3m3", "RH2M3m3")],
  ra.ls[["ra1m"]][ , c("PRECTOTCORR1m", "T2M1m", "RH2M1m",
                 "PRECTOTCORR3m", "T2M3m", "RH2M3m")]
))
# confirmed perfect correlations between two datasets for common samples
round(cor(cbind(
  ra3mSub1m[ , c("PRECTOTCORR1m3", "T2M1m3", "RH2M1m3",
                 "PRECTOTCORR3m3", "T2M3m3", "RH2M3m3")],
  ra.ls[["ra1m"]][ , c("PRECTOTCORR1m", "T2M1m", "RH2M1m",
                       "PRECTOTCORR3m", "T2M3m", "RH2M3m")]
), method = "kendall"), 2)

# output files with relative abundance response and covariates ------------------

saveRDS(ra.ls[["ra3m"]], file = "ra3m_obs_covs_stra.rds")
saveRDS(ra.ls[["ra1m"]], file = "ra1m_obs_covs_stra.rds")

