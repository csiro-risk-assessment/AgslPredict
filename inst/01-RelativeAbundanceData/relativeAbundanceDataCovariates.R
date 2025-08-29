# VectorBase relative abundance data and covariates ----------------------------
# zipping together VectorBase relative abundance data
# compiled to grid cell resolution

# read in relative abundance data previously output to this directory
# by running "RelativeAbundanceDatePrep.R" also in this directory
ra <- readRDS("VB_PCR_rel_abundance_stra.rds")

# non meteorological covariates ------------------------------------------------

# Lambert azimuthal equal area
# see Steinwand et al. (1995)

library(terra)
active.r <- rast("../covariates_spatial/activeAfrica.tif")
a.crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
same.crs(crs(active.r), a.crs) # TRUE
a.crs <- crs(active.r)

## read in elevation -----------------------------------------------------------

elev.r <- rast("../covariates_spatial/elev.tif")

## read in remaining spatial covariate files ------------------------------

d2c.r <- rast("../covariates_spatial/d2c.tif")
ppw.r <- rast("../covariates_spatial/lakesrivers.tif")

res(d2c.r); res(ppw.r); res(elev.r)

## assign spatial covariate values to relative abundance data ------------------

# project coords from abundance data and find corresponding cells
# (could also use cell id but sticking with geographic overlay)
coords <- cbind(ra$lon.centroid, ra$lat.centroid)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, a.crs)
par(mfrow = c(1, 1))
plot(active.r)
plot(proj.coords, add = TRUE)

o.d2c <- extract(d2c.r, proj.coords)
sum(o.d2c[ , 2] >= 0) == nrow(proj.coords)
o.ppw <- extract(ppw.r, proj.coords)
sum(o.ppw[ , 2] >= 0) == nrow(proj.coords)
o.elev <- extract(elev.r, proj.coords) 
sum(o.elev[ , 2] >= -155) == nrow(proj.coords)

ra$d2c <- o.d2c[ , 2]
ra$ppw <- o.ppw[ , 2]
ra$elev <- o.elev[ , 2]

#  vector intervention covariates (spatio-temporal) ----------------------------

Y <- 2001:2021

## IRS -------------------------------------------------------------------------

ra$IRS <- NA
for (y in Y) {
  irs <- rast(paste0("../vector_intervention/IRS_", y, ".tif"))
  same.crs(active.r, irs) # TRUE
  o.irs.y <- extract(irs, proj.coords)
  ra$IRS[ra$origin.year == y] <- o.irs.y[ra$origin.year == y, 2]
}

### IRS misalignment near coast ------------------------------------------------

# IRS intervention data missing from these four cells near coast
na.irs <- ra$cell[which(is.na(ra$IRS))] # "1133068", "1133068", "1256195", "644773"

# missing data occur in four different years from 2001 to 2004:
ra[ra$cell%in%na.irs, ]

# for these misaligned cells assign nearest neighbour intervention data in space and time
# surrounding regions show zero intervention:
pt <- crds(proj.coords[ra$cell == 1256195, ])
irs <- rast(paste0("../vector_intervention/IRS_", 2001, ".tif"))
plot(irs, xlim = c(-1e6, 1e6) + pt[1, "x"], ylim = c(-1e6, 1e6) + pt[1, "y"])
plot(proj.coords[ra$cell == 1256195, ], add = TRUE, col = 'red') # zero
ra[ra$cell == 1256195, "IRS"] <- 0

pt <- crds(proj.coords[ra$cell == 1133068, ])
irs <- rast(paste0("../vector_intervention/IRS_", 2002, ".tif"))
plot(irs, xlim = c(-1e6, 1e6) + pt[1, "x"], ylim = c(-1e6, 1e6) + pt[1, "y"])
plot(proj.coords[ra$cell == 1256195, ], add = TRUE, col = 'red') # zero
irs <- rast(paste0("../vector_intervention/IRS_", 2003, ".tif"))
plot(irs, xlim = c(-1e6, 1e6) + pt[1, "x"], ylim = c(-1e6, 1e6) + pt[1, "y"])
plot(proj.coords[ra$cell == 1256195, ], add = TRUE, col = 'red') # zero
ra[ra$cell == 1133068, "IRS"] <- 0

pt <- crds(proj.coords[ra$cell == 644773, ])
irs <- rast(paste0("../vector_intervention/IRS_", 2004, ".tif"))
plot(irs, xlim = c(-1e6, 1e6) + pt[1, "x"], ylim = c(-1e6, 1e6) + pt[1, "y"])
plot(proj.coords[ra$cell == 644773, ], add = TRUE, col = 'red') # zero
ra[ra$cell == 644773, "IRS"] <- 0

## ITN -------------------------------------------------------------------------

ra$ITN <- NA
for (y in Y) {
  itn <- rast(paste0("../vector_intervention/ITN_", y, ".tif"))
  same.crs(active.r, itn) # TRUE
  o.itn.y <- extract(itn, proj.coords)
  ra$ITN[ra$origin.year == y] <- o.itn.y[ra$origin.year == y, 2]
}

# human population density -----------------------------------------------------

ra$hpd <- NA
for (y in Y) {
  hpd <- rast(paste0("../human_pop/human_density_", y, ".tif"))
  same.crs(active.r, hpd) # TRUE
  o.hpd.y <- extract(hpd, proj.coords)
  ra$hpd[ra$origin.year == y] <- o.hpd.y[ra$origin.year == y, 2]
}

# human pop density data missing
na.hpd <- ra$cell[which(is.na(ra$hpd))] # 1402901 in year 2001

pt <- crds(proj.coords[ra$cell == 1402901, ])
hpd <- rast(paste0("../human_pop/human_density_", 2001, ".tif"))
plot(hpd, xlim = c(-1e5, 1e5) + pt[1, "x"], ylim = c(-1e5, 1e5) + pt[1, "y"])
plot(proj.coords[ra$cell == 1402901, ], add = TRUE, col = 'red') # coastal ocean
# assign nearest coastal density
points(3203366 - 5000*2.75, -2414849 + 3500, col = 'cyan')
extract(hpd, matrix(c(3203366 - 5000*2.75, -2414849 + 3500), nrow = 1), cells = TRUE)
ra[ra$cell == 1402901, "hpd"] <- 11.69416

# meteorological covariates ----------------------------------------------------

# test model: two datasets
# 31 day max duration for sample (~1 month)

ra1m <- ra

ra.names <- c("ra1m")
ra.ls <- list(ra1m = NULL)

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
hist(ra.ls[["ra1m"]]$PRECTOTCORR1m)
hist(ra.ls[["ra1m"]]$T2M1m)
hist(ra.ls[["ra1m"]]$RH2M1m)
hist(ra.ls[["ra1m"]]$PRECTOTCORR3m)
hist(ra.ls[["ra1m"]]$T2M3m)
hist(ra.ls[["ra1m"]]$RH2M3m)

# vector intervention summary stats --------------------------------------------

# Malaria Atlas Project:
# "Proportion of households covered with Indoor Residual Spraying during a defined year 2000-2022"
tIRS <- c(tapply(ra.ls[["ra1m"]]$IRS, list(ra.ls[["ra1m"]]$cell), "max", na.rm = TRUE))
sum(tIRS > 0.1)/length(tIRS)

# Malaria Atlas Project:
# "Proportion of population that sleeps under an Insecticide-Treated Net during a defined year 2000-2022"
tITN <- c(tapply(ra.ls[["ra1m"]]$ITN, list(ra.ls[["ra1m"]]$cell), "max", na.rm = TRUE))
sum(tITN > 0.1)/length(tITN)

# output files with relative abundance response and covariates ------------------

saveRDS(ra.ls[["ra1m"]], file = "ra1m_obs_covs_stra.rds")

