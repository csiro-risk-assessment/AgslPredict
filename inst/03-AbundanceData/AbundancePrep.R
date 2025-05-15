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

# read in spatial scope of RA 
library(terra)
active.r <- rast("../covariates_spatial/activeAfrica.tif")
a.crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" 
same.crs(crs(active.r), a.crs) # TRUE
a.crs <- crs(active.r)

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
proj.coords <- project(v, a.crs)

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
or.crds <- vect(or.crds, crs = a.crs)
por.crds <- project(or.crds, "epsg:4326")
o.cells <- cbind(o.cells, crds(por.crds))

a.dat$cell <- o.cells[ , "cell"]
a.dat$lon.centroid <- o.cells[ , "x"]
a.dat$lat.centroid <- o.cells[ , "y"]

plot(a.dat$Longitudes, a.dat$Latitudes)
points(a.dat$lon.centroid, a.dat$lat.centroid, pch = 3, col = 'red')

# there is a single study pre-2001
unique(a.dat[lubridate::year(lubridate::date(a.dat$Collection.date.range)) < 2001, "cell"])
# VBA2882243
# DOI:10.1186/1471-2334-12-275,DOI:10.1186/1756-3305-6-10
# cell: 684641
points(a.dat$lon.centroid[a.dat$cell == 684641], 
       a.dat$lat.centroid[a.dat$cell == 684641], 
       pch = 19, col = 'blue')

# this is city of Mbandjock with ITN intervention: Antonio-Nkondjio et al. (2013)
unique(a.dat[lubridate::year(lubridate::date(a.dat$Collection.date.range)) < 2001, c("Longitudes", "Latitudes")])
4 + 27/60 # lat Mbandjock, Cameroon
11 + 54/60 # long Mbandjock, Cameroon

# collection methods for analysis ----------------------------------------------

# collection method management
# exclude larvae collections - targeting adult females
a.dat <- a.dat[a.dat$Collection.protocols != "collection of larvae from dippers", ]

# drop all records with collection methods having less than 10 occurrences
rareMethods <- names(table(a.dat$Collection.protocols))[table(a.dat$Collection.protocols) <= 10]

# collection method management
# exclude larvae collections - targeting adult females
dat <- dat[dat$Collection.protocols != "collection of larvae from dippers", ]
# drop all records with collection methods having less than 10 occurrences
# defined as date by cell by protocol combinations
rareMethods <- c(
  "animal biting catch - indoors", "catch of live specimens",
  "house resting - daytime catch", "swarm net catch"
)
a.dat <- a.dat[!a.dat$Collection.protocols%in%rareMethods, ]

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

all.equal(sort(unique(a.dat[, "Projects"])), unique(sort(Abund.cit[ , "Projects"])))

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
proj.coords <- project(v, a.crs)

plot(active.r)
plot(proj.coords, add = TRUE, col = "darkorange")


