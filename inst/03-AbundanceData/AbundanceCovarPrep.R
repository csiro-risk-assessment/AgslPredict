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

library(terra)
active.r <- rast("../covariates_spatial/activeAfrica.tif")
a.crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
same.crs(crs(active.r), a.crs) # TRUE
a.crs <- crs(active.r)

res(active.r) # 5000 5000

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

d2c.r <- rast(d2c[nrow(d2c):1, ], crs = a.crs,
              extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
d2r.r <- rast(d2r[nrow(d2r):1, ], crs = a.crs,
              extent =  c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
pop.r <- rast(pop[nrow(pop):1, ], crs = a.crs,
              extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))

res(d2c.r); res(d2r.r); res(pop.r); res(elev.r)

# assign covariate values to abundance data ------------------------------------

# project coords from abundance data and find corresponding cells
# (could also use cell id but sticking with geographic overlay)
coords <- cbind(abund$lon.centroid, abund$lat.centroid)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, a.crs)

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

## met data covariates ---------------------------------------------------------

abundance <- cbind(abund, pabund[ , c(names(lag.p), names(lag.t), names(lag.r))])

## vector intervention ---------------------------------------------------------

# years with observations
Y <- sort(unique(lubridate::year(as.Date(abundance$date))))
Y.pre2001 <- Y[Y < 2001]
Y <- Y[Y >= 2001]

abundance$year <- lubridate::year(as.Date(abundance$date))

## IRS -------------------------------------------------------------------------

abundance$IRS <- NA
for (y in Y) {
  irs <- rast(paste0("../vector_intervention/IRS_", y, ".tif"))
  same.crs(active.r, irs) # TRUE
  o.irs.y <- extract(irs, proj.coords)
  abundance$IRS[abundance$year == y] <- o.irs.y[abundance$year == y, 2]
}

### IRS missing for early years-------------------------------------------------

# this is city of Mbandjock with ITN intervention: Antonio-Nkondjio et al. (2013)
# described in DOI:10.1186/1756-3305-6-10 as having
# no IRS treatment since 1960s
unique(abundance[is.na(abundance$IRS), "cell"]) # 684641
unique(abundance[is.na(abundance$IRS), "year"]) # 1997 1998

abundance$IRS[is.na(abundance$IRS)] <- 0

## ITN -------------------------------------------------------------------------

abundance$ITN <- NA
for (y in Y) {
  itn <- rast(paste0("../vector_intervention/ITN_", y, ".tif"))
  same.crs(active.r, itn) # TRUE
  o.itn.y <- extract(itn, proj.coords)
  abundance$ITN[abundance$year == y] <- o.itn.y[abundance$year == y, 2]
}

### ITN missing for early years ------------------------------------------------

# this is city of Mbandjock with ITN intervention: Antonio-Nkondjio et al. (2013)
# described in DOI:10.1186/1756-3305-6-10
# "January to June 1997 represented the period before bed net coverage and September 1997 to
# September 1998 was the period after bed net coverage"
# Table 1:
# 3400 total beds in Mbandjock
# 31 pre-existing nets
# 2454 nets distributed
unique(abundance[is.na(abundance$ITN), "cell"]) # 684641
unique(abundance[is.na(abundance$ITN), "year"]) # 1997 1998

unique(abundance[abundance$cell == 684641 & abundance$date < "1997-07-01", "date"])
unique(abundance[abundance$cell == 684641 & abundance$date >= "1997-07-01", "date"])

abundance[abundance$cell == 684641 & abundance$date < "1997-07-01", "ITN"] <- 31/3400
abundance[abundance$cell == 684641 & abundance$date >= "1997-07-01", "ITN"] <- 2454/3400

# vector intervention summary stats --------------------------------------------

# Malaria Atlas Project:
# "Proportion of households covered with Indoor Residual Spraying during a defined year 2000-2022"
tIRS <- c(tapply(abundance$IRS, list(abundance$cell), "max", na.rm = TRUE))
sum(tIRS > 0.1)/length(tIRS)

# Malaria Atlas Project:
# "Proportion of population that sleeps under an Insecticide-Treated Net during a defined year 2000-2022"
tITN <- c(tapply(abundance$ITN, list(abundance$cell), "max", na.rm = TRUE))
sum(tITN > 0.1)/length(tITN)

# save abundance file ----------------------------------------------------------

saveRDS(abundance, file = "abund_obs_covs.rds")
