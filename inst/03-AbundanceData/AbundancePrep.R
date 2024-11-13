# download abundance data from VectorBase
# and prepare for analysis

# download abundance data ------------------------------------------------------

# increase timeout for internet download
options(timeout = max(1200, getOption("timeout")))

download.file("https://vectorbase.org/common/downloads/popbio-map-legacy/VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv.gz",
              destfile = "VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv.gz")

# check download
file.exists("VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv.gz")
# unzip using R.utils
if(!(file.exists("VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv"))) {
  R.utils::gunzip("VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv.gz", remove = FALSE)}
# read in data using data.table
dat <- data.table::fread("VectorBase-popbio-map-rel64-July2023-Abundance-with-zeroes.csv",
                         header = TRUE, stringsAsFactors = FALSE)
dat <- as.data.frame(dat)

col.names <- colnames(dat)
col.names <- gsub(" ", ".", col.names)
col.names <- gsub("\\(", ".", col.names)
colnames(dat) <- gsub("\\)", ".", col.names)

# subset abundance data to bounding box ----------------------------------------

# read in spatial scope of RA in projected CRS
active <- read.csv("../covariates_spatial/active.csv", skip = 1, header = FALSE) # or skip = 2?
#active <- read.csv("inst/covariates_spatial/active.csv", skip = 1, header = FALSE)
all.equal(dim(active), c(1280, 1520)) # expecting same number of rows and columns for all covariate data
active <- as.matrix(active)

library(terra)
# Lambert azimuthal equal area
# see Steinwand et al. (1995)
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

plot(active.r, ext = c(-160000 - 20000*2, -140000 + 20000*2, -69000 - 40000*2, -68000 + 2*40000))

project(ext(active.r), crs(active.r), "epsg:4326")

dat <- dat[dat$Longitudes > -23 & dat$Longitudes < 57, ]
dat <- dat[dat$Latitudes > -34 & dat$Latitudes < 25, ]

# data management -------------------------------------------------------------

# Project and study https://doi.org/10.1186/1475-2875-9-187 is
# indoor CDC light trap
table(dat$Collection.protocols[dat$Projects == "VBP0000614"])
unique(dat$Collection.protocols)
# relabel as "indoor light trap catch"
dat[dat$Projects == "VBP0000614", "Collection.protocols"] <- "indoor light trap catch"

# small proportion of collections report "catch of live specimens" instead of
# "indoor light trap catch" for project VBP0000774
# no citation info
# drop catch of live specimen records for this project
table(dat[dat$Project == "VBP0000774", "Collection.protocols"])
unique(dat[dat$Project == "VBP0000774", "Citations"])
dat <- dat[!(dat$Project == "VBP0000774" & dat$Collection.protocols == "catch of live specimens"), ]
# project VBP0000624  also reports "catch of live specimens"
table(dat[dat$Collection.protocols == "catch of live specimens", "Projects"])
# source "DOI:10.1186/1475-2875-13-27" reports combo of outdoor light traps and hlc
# also larval collection but those samples don't include adults. Retain.

# species ids assigned to An. gambiae s.l.
agsl.id <- c(
  "Anopheles arabiensis",
  "Anopheles coluzzii",
  "Anopheles gambiae sensu lato",
  "Anopheles gambiae sensu stricto",
  "Anopheles gambiae x Anopheles coluzzii"
)
dat <- dat[dat$Species%in%agsl.id, ]

# about 4% of biting catch is scored male, 95% female and remaining mixed or unknown
table(dat$Sex[grepl("biting catch", dat$Collection.protocols)])/sum(grepl("biting catch", dat$Collection.protocols))

# subset to adult females
dat <- dat[dat$Sex == "female" & dat$Developmental.stage == "adult", ]

# what is a sample? It is a unique "Collection.ID"
# VB data dictionary:
# "Collection ID Unique ID for all organisms collected at a single time. For example, all the
# organisms collected in a single passive trap during one collection period would share the
# same Collection ID"
dat$samp.id <- paste(dat$Latitudes, dat$Longitudes, dat$Collection.date.range,
                     dat$Collection.protocols, dat$Projects)
samps <- sort(unique(dat$samp.id))

# create dataframe of species frequencies
# roll through each sample
# preserve metadata from first record that appears in each sample
for (i in 1:length(samps)) {

  temp.a <- dat[dat$samp.id == samps[i], ]
  temp.a$units <- length(unique(temp.a$Collection.ID))
  temp.a$Agsl <- sum(temp.a$Specimens.collected)
  temp.a$Species <- temp.a$Sample.ID <- temp.a$Label <- NULL # drop sample specific info
  temp.a <- temp.a[1, ]

  if (i == 1) {
    # initialise abundance data set
    a.dat <- temp.a[1, ]
  } else {
    a.dat <- rbind(a.dat, temp.a)
  }
  cat(i, " out of ", length(samps), " completed ", "\n")
}

# 2% of observations have indeterminate date of sampling
sum(table(nchar(a.dat$Collection.date.range))[c("7", "18")])/nrow(a.dat)
# subset to observations with defined dates
a.dat <- a.dat[nchar(a.dat$Collection.date.range) == 10 | nchar(a.dat$Collection.date.range) == 21, ]

# assign overnight catches, e.g, HLC, to first day of sampling
split.date <- strsplit(a.dat$Collection.date.range, "/")
a.dat$date <- sapply(split.date, function(x) x[1])

# visualise occurrence locations
plot(a.dat$Longitudes[a.dat$Agsl > 0], a.dat$Latitudes[a.dat$Agsl > 0],
     xlim = range(a.dat$Longitudes), ylim = range(a.dat$Latitudes),
     main = "VectorBase abundance",
     xlab = "long", ylab = "lat")

# aggregate to grid-------------------------------------------------------------

# project coords from abundance data and find corresponding cells
coords <- cbind(a.dat$Longitudes, a.dat$Latitudes)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, crs)

plot(active.r)
plot(proj.coords, add = TRUE, col = 'red')
oa <- extract(active.r, proj.coords)
# not all observation coords correspond to an active cell
sum(oa[ , 2] == 1) == nrow(proj.coords)
plot(proj.coords[oa[ , 2] == 0, ], add = TRUE, col = 'green', pch = 3)

unique(paste(crds(proj.coords)[oa[ , 2] == 0, 1], crds(proj.coords)[oa[ , 2] == 0, 2]))
plot(proj.coords[oa[ , 2] == 0, ])

o.cells <- cells(active.r, proj.coords)
# convert centroid of overlaid cells to lat/long
r.crds <- crds(active.r)
or.crds <- r.crds[o.cells[ , "cell"], ]
# project to lat/long
or.crds <- vect(or.crds, crs = crs)
por.crds <- project(or.crds, "epsg:4326")
o.cells <- cbind(o.cells, crds(por.crds))

a.dat$cell <- o.cells[ , "cell"]
a.dat$lon.centroid <- o.cells[ , "x"]
a.dat$lat.centroid <- o.cells[ , "y"]

plot(a.dat$Longitudes, a.dat$Latitudes)
points(a.dat$lon.centroid, a.dat$lat.centroid, pch = 3, col = 'red')

# abundance data source map ----------------------------------------------------
library(sf)
library(rnaturalearth)

africa <- ne_countries(continent = 'Africa', scale = 'medium', returnclass = 'sf')
same.crs(africa, proj.coords)

proj.coords.newcrs <- st_transform(st_as_sf(proj.coords), st_crs(africa))
same.crs(africa, proj.coords.newcrs)

library(ggplot2)
p <- ggplot(data = africa) +
  geom_sf() +
  theme_bw() +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  ggtitle("Location of abundance data sets") +
  geom_sf(
    data = proj.coords.newcrs,
    size = 1,
    shape = 16,
    colour = "blue"
  )
p
ggsave(
  "AbundDataMap.png",
  p,
  width = 20,
  height = 20,
  units = "cm"
)

# sample compilation -----------------------------------------------------------

# given aggregation to cells
# what is a sample? It is a unique "Collection.ID"
# VB data dictionary:
# "Collection ID Unique ID for all organisms collected at a single time. For example, all the
# organisms collected in a single passive trap during one collection period would share the
# same Collection ID"
# sample: same cell id, same date, same collection protocol. Project agnostic.
a.dat$samp.id <- paste(a.dat$cell, a.dat$date, a.dat$Collection.protocols)
samps <- sort(unique(a.dat$samp.id))

# create dataframe of species frequencies
# roll through each sample
# preserve metadata from first record that appears in each sample
for (i in 1:length(samps)) {

  temp.a <- a.dat[a.dat$samp.id == samps[i], ]
  temp.a$units <- sum(temp.a$units)
  temp.a$Agsl <- sum(temp.a$Agsl)
  temp.a <- temp.a[1, ]

  if (i == 1) {
    # initialise abundance data set
    abundance <- temp.a[1, ]
  } else {
    abundance <- rbind(abundance, temp.a)
  }

}

## Create a table of the citations for the abundance data ----------------------
Abund.cit <- unique(abundance[, c("Projects", "Citations")])

## Export as a table
library(xtable)
print(
  xtable(Abund.cit, caption = "Caption here"),
  include.rownames = F,
  size = "small",
  file = "AbundCit.txt"
)

# clean up abundance -----------------------------------------------------------
abundance$Record.type <- abundance$Sample.type <- abundance$Collection.ID <-
  abundance$Collection.date.range <- abundance$Projects <-
  abundance$Latitudes <- abundance$Longitudes <- abundance$Locations <-
  abundance$Citations <- abundance$Tag <- abundance$Attractants <-
  abundance$Usage.license <- abundance$Geolocation.provenance <-
  abundance$Geolocation.precision <- NULL

# Export abundance data
saveRDS(abundance, file = "abundance.rds")

# visualise spatial
coords <- cbind(abundance$lon.centroid, abundance$lat.centroid)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, crs)

plot(active.r)
plot(proj.coords, add = TRUE, col = "darkorange")



# # abundance aggregated to grid cells
# abund <- readRDS("abundance.rds")

# # weekly running means of precip -----------------------------------------------
# # lagged from one day to duration of immature life history stage
# # for each observation
#
# # lat/long for each cell
# precip.cells <- unique(abund$cell)
# precip.lon <- precip.lat <- numeric(length(precip.cells))
# for (i in 1:length(precip.cells)) {
#   precip.lon[i] <- unique(abund[abund$cell == precip.cells[i], "lon.centroid"])
#   precip.lat[i] <- unique(abund[abund$cell == precip.cells[i], "lat.centroid"])
# }
# precip.locations <- as.data.frame(cbind(precip.cells, precip.lon, precip.lat))
#
# imm.duration <- 10 # duration of immature phase
# p.duration <- 7 # duration of precip running mean
#
# # running weekly averages
# # for lag of one day to lag of imm.duration days wrt every observation
# lag.p <- lag.t <- lag.r <- as.data.frame(matrix(NA, nrow = nrow(abund), ncol = imm.duration,
#                                                 dimnames = list(NULL, paste0("lag", 1:imm.duration))
# ))
# names(lag.p) <- paste0("p_", names(lag.p))
# names(lag.t) <- paste0("t_", names(lag.t))
# names(lag.r) <- paste0("r_", names(lag.r))
#
# pabund <- cbind(abund, lag.p, lag.t, lag.r)
#
# ## precipitations
# for (i in 1:nrow(precip.locations)) {
#
#   rec <- read.csv(paste0("../meteo_for_abundance/precip_", precip.cells[i], ".csv"),
#                   header = TRUE)
#
#   # subset dataframe to correspond with cell
#   sub.cell <- pabund[pabund$cell == precip.cells[i], ]
#
#   # lagged weekly averages of precip
#   for (j in 1:nrow(sub.cell)) {
#     id <- which(rec$YYYYMMDD == sub.cell$date[j])
#     # lags from one day to imm.duration days
#     for (k in 1:imm.duration) {
#       # p.duration running means over lags of one day to imm.duration days
#       sub.cell[j, colnames(lag.p)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "PRECTOTCORR"])
#       sub.cell[j, colnames(lag.t)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "T2M"])
#       sub.cell[j, colnames(lag.r)[k]] <- mean(rec[(id - k - p.duration + 1):(id - k), "RH2M"])
#     }
#   }
#
#   pabund[pabund$cell == precip.cells[i], ] <- sub.cell
#
# }
#
# all.equal(abund, pabund[ , names(abund)])
#
# # spatial covariates -----------------------------------------------------------
#
# # following uses same code as "relativeAbundancePrep.R"
# # switching data.frame ra for abund
#
# # Lambert azimuthal equal area
# # see Steinwand et al. (1995)
#
# crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
#
# # use  limits denoted in first line of csv files
# # should read: #xmin=-4099134.0	ymin=-4202349.0	cell_size=5000.0	nx=1520	ny=1280
#
# # check active cells for correct extent info
# all(read.csv("../covariates_spatial/active.csv", header = FALSE, nrow = 1) ==
#       c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# # read in spatial scope of RA in projected CRS
# active <- read.csv("../covariates_spatial/active.csv", skip = 1, header = FALSE)
# all.equal(dim(active), c(1280, 1520)) # expecting same number of rows and columns for all covariate data
# active <- as.matrix(active)
# library(terra)
# # extent and crs
# active.r <- rast(active[nrow(active):1, ], crs = crs,
#                  extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
# res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above
#
# # read in elevation
# elev.r <- rast("../covariates_spatial/elev.tif")
#
# # read in covariates with existing .csv files (10 August 2023)
#
# # check for correct extent info
# all(read.csv("../covariates_spatial/d2c.csv", header = FALSE, nrow = 1) ==
#       c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# all(read.csv("../covariates_spatial/d2r.csv", header = FALSE, nrow = 1) ==
#       c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
# all(read.csv("../covariates_spatial/pop.csv", header = FALSE, nrow = 1) ==
#       c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
#
# d2c <- read.csv("../covariates_spatial/d2c.csv", skip = 1, header = FALSE)
# d2r <- read.csv("../covariates_spatial/d2r.csv", skip = 1, header = FALSE)
# pop <- read.csv("../covariates_spatial/pop.csv", skip = 1, header = FALSE)
#
# all.equal(dim(d2c), c(1280, 1520))
# all.equal(dim(d2r), c(1280, 1520))
# all.equal(dim(pop), c(1280, 1520))
#
# d2c <- as.matrix(d2c)
# d2r <- as.matrix(d2r)
# pop <- as.matrix(pop)
#
# d2c.r <- rast(d2c[nrow(d2c):1, ], crs = crs,
#               extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
# d2r.r <- rast(d2r[nrow(d2r):1, ], crs = crs,
#               extent =  c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
# pop.r <- rast(pop[nrow(pop):1, ], crs = crs,
#               extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
#
# res(d2c.r); res(d2r.r); res(pop.r); res(elev.r)
#
# # assign covariate values to abundance data ------------------------------------
#
# # project coords from abundance data and find corresponding cells
# # (could also use cell id but sticking with geographic overlay)
# coords <- cbind(abund$lon.centroid, abund$lat.centroid)
# v <- vect(coords, crs = "epsg:4326")
# proj.coords <- project(v, crs)
#
# plot(active.r)
# plot(proj.coords, add = TRUE)
#
# o.d2c <- extract(d2c.r, proj.coords)
# o.d2r <- extract(d2r.r, proj.coords)
# o.pop <- extract(pop.r, proj.coords) # no missing values
# o.elev <- extract(elev.r, proj.coords) # no missing values
#
# abund$d2c <- o.d2c[ , 2]
# abund$d2r <- o.d2r[ , 2]
# abund$pop <- o.pop[ , 2]
# abund$elev <- o.elev[ , 2]
#
# # bind spatial and temporal covariates -----------------------------------------
#
# abundance <- cbind(abund, pabund[ , c(names(lag.p), names(lag.t), names(lag.r))])
# saveRDS(abundance, file = "abund_obs_covs.rds")
