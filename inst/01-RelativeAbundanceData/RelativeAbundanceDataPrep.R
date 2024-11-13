library(dplyr)

# download and import data -----------------------------------------------------
#raw.dir <- "./"
download.file("https://vectorbase.org/common/downloads/popbio-map-legacy/VectorBase-popbio-map-rel64-July2023-Samples.csv.gz",
              destfile = "VectorBase-popbio-map-rel64-July2023-Samples.csv.gz")
# check download
file.exists("VectorBase-popbio-map-rel64-July2023-Samples.csv.gz")
# unzip using R.utils
if(!(file.exists("VectorBase-popbio-map-rel64-July2023-Samples.csv"))) {
R.utils::gunzip("VectorBase-popbio-map-rel64-July2023-Samples.csv.gz", remove = FALSE)}
# read in data using data.table
dat <- data.table::fread("VectorBase-popbio-map-rel64-July2023-Samples.csv",
                         header = TRUE, stringsAsFactors = FALSE)
dat <- as.data.frame(dat)
col.names <- colnames(dat)
colnames(dat) <- gsub(" ", ".", col.names)

# species identification -------------------------------------------------------

# account for chormosomal forms
dat$sp <- dat$Species
dat$sp <- ifelse(dat$sp == "Anopheles gambiae chromosomal form Bamako",
                 "Anopheles gambiae sensu stricto",
                 dat$sp)
dat$sp <- ifelse(dat$sp == "Anopheles gambiae chromosomal form Savanna",
                 "Anopheles gambiae sensu stricto",
                 dat$sp)
dat$sp <- ifelse(dat$sp == "Anopheles gambiae chromosomal form Mopti",
                 "Anopheles coluzzii",
                 dat$sp)

# subset to datasets that identify at least one of following species in "sp"
# a single project VBP0000812 reports "Anopheles gambiae x Anopheles coluzzii"
# drop from species composition analysis
sp.id <- c("Anopheles arabiensis",
           "Anopheles coluzzii",
           "Anopheles gambiae sensu stricto")
dat <- dat[dat$sp%in%sp.id, ]

# geographical and temporal filtering ------------------------------------------

# make sure we are only taking data in Africa / Arabian peninsula
dat <- dat[dat$Latitudes > -34, ]
dat <- dat[dat$Latitudes < 20, ]
dat <- dat[dat$Longitudes > -20, ]
dat <- dat[dat$Longitudes < 52, ]


# remove records with no year or no month information or duration greater than 3 months
dat <- dat[nchar(dat$Collection.date.range) != 0, ] # no dates
dat <- dat[nchar(dat$Collection.date.range) != 4, ] # only year
dat <- dat[nchar(dat$Collection.date.range) != 9, ] # only year range
dat <- dat[nchar(dat$Collection.date.range) != 14, ] # only year collection
dat <- dat[nchar(dat$Collection.date.range) != 23, ] # duration beyond three months for these records
dat <- dat[nchar(dat$Collection.date.range) != 31, ] # duration beyond three months for these records
dat <- dat[nchar(dat$Collection.date.range) != 47, ] # duration beyond three months for these records
dat <- dat[nchar(dat$Collection.date.range) != 63, ] # duration beyond three months for these records

library(zoo)
library(lubridate)

dat$date1 <- dat$date2 <- NA
for (d in unique(dat$Collection.date.range)) {

  if (nchar(d) == 7) {
    # recorded to month
    month <- as.yearmon(d)
    dat$date1[dat$Collection.date.range == d] <- paste0(d, "-01")
    dat$date2[dat$Collection.date.range == d] <- paste0(d, "-", days_in_month(month))
  }
  if (nchar(d) == 10) {
    # recorded to day
    dat$date1[dat$Collection.date.range == d] <-
      dat$date2[dat$Collection.date.range == d] <- d
  }
  if (nchar(d) == 15) {
    # duration between months
    splits <- strsplit(d, "/")[[1]]
    months <- as.yearmon(splits)
    dat$date1[dat$Collection.date.range == d] <- paste0(splits[1], "-01")
    dat$date2[dat$Collection.date.range == d] <- paste0(splits[2], "-", days_in_month(months[2]))
  }
  if (nchar(d) == 18) {
    # duration between month and spurious date equal to first of same month
    splits <- strsplit(d, "/")[[1]]
    month <- as.yearmon(splits[1])
    dat$date1[dat$Collection.date.range == d] <- paste0(splits[1], "-01")
    dat$date2[dat$Collection.date.range == d] <- paste0(splits[1], "-", days_in_month(month))
  }
  if (nchar(d) == 21) {
    # duration between two dates
    splits <- strsplit(d, "/")[[1]]
    dat$date1[dat$Collection.date.range == d] <- splits[1]
    dat$date2[dat$Collection.date.range == d] <- splits[2]
  }

}

dat$diffDate <- as.numeric(date(dat$date2) - date(dat$date1))

# percentage reported to month
sum(nchar(dat$Collection.date.range) == 7)/nrow(dat)
# percentage reported to day
sum(nchar(dat$Collection.date.range) == 10)/nrow(dat)
# percentage reported to date range
sum(nchar(dat$Collection.date.range) == 21)/nrow(dat)
# of those collections with duration between dates
sum(dat$diffDate[nchar(dat$Collection.date.range) == 21] <= 1)/nrow(dat[nchar(dat$Collection.date.range) == 21, ])
# samples with diffDate greater than 31 days
sum(dat$diffDate > 31)/nrow(dat)
# samples with diffDate greater than 93 days (about 3 months)
sum(dat$diffDate > 93)/nrow(dat)

# define year of origin for sample
dat$origin.year <- year(dat$date1)
# remove all dates with samples initiated before 2001
dat <- dat[dat$origin.year >= 2001, ]

# include samples with diffDate not greater than ~3 months or 93 days (allowing for overnight collections)
dat <- dat[dat$diffDate <= 93, ]

# species frequencies ----------------------------------------------------------

# defined sample id
# identified by shared lat/long, date collection protocol, and available data types
# note: lots of examples of mixed collections sources, e.g., both indoors and outdoors
# missing collection.data.range records omitted above
dat$samp.id <- paste(dat$Latitudes, dat$Longitudes, dat$Collection.date.range,
                     dat$Collection.protocols, dat$Projects)
samps <- sort(unique(dat$samp.id))

# create dataframe of species frequencies
# roll through each sample
# preserve metadata from first record that appears in each sample
for (i in 1:length(samps)) {

  temp.ra <- dat[dat$samp.id == samps[i], ]
  temp.ra$Aa <- temp.ra$Ac <- temp.ra$Ag <- NA
  temp.ra$Aa <- sum(temp.ra$sp == "Anopheles arabiensis")
  temp.ra$Ac <- sum(temp.ra$sp == "Anopheles coluzzii")
  temp.ra$Ag <- sum(temp.ra$sp == "Anopheles gambiae sensu stricto")
  temp.ra$Species <- temp.ra$sp <- NULL # drop old species info
  temp.ra <- temp.ra[1, ]

  if (i == 1) {
    # initialise relative abundance data set
    ra.dat <- temp.ra[1, ]
  } else {
    ra.dat <- rbind(ra.dat, temp.ra)
  }

}

# visualise occurrence locations by species
plot(ra.dat$Longitudes[ra.dat$Aa > 0], ra.dat$Latitudes[ra.dat$Aa > 0],
     xlim = range(ra.dat$Longitudes), ylim = range(ra.dat$Latitudes),
     main = "VectorBase: Presences by species (PCR protocol Legacy download)",
     xlab = "long", ylab = "lat")
points(ra.dat$Longitudes[ra.dat$Ag > 0], ra.dat$Latitudes[ra.dat$Ag > 0], pch = 4, col = 'green')
points(ra.dat$Longitudes[ra.dat$Ac > 0], ra.dat$Latitudes[ra.dat$Ac > 0], pch = 3, col = 'blue')
legend("bottomleft", pch = c(1, 3, 4), col = c('black', 'blue', "green"),
       legend = c("An. arabiensis", "An. coluzzii", "An. gambiae sensu stricto"))

# collation to grid ----------------------------------------------------------

# collate to grid
active <- read.csv("../covariates_spatial/active.csv", skip = 1, header = FALSE)
all.equal(dim(active), c(1280, 1520)) # expected number of rows and columns
active <- as.matrix(active)
library(terra)
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
# extent and crs from active.csv
active.r <- rast(active[nrow(active):1, ], crs = crs,
                 extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))
res(active.r) # 5000 5000 with skip = 1 but 5000.000 5003.909 with skip = 2 in read.csv above

# project coords from abundance data and find corresponding cells
coords <- cbind(ra.dat$Longitudes, ra.dat$Latitudes)
v <- vect(coords, crs = "epsg:4326")
proj.coords <- project(v, crs)
oa <- extract(active.r, v)

o.cells <- cells(active.r, proj.coords)
# convert centroid of overlaid cells to lat/long
r.crds <- crds(active.r)
or.crds <- r.crds[o.cells[ , "cell"], ]
# project to lat/long
or.crds <- vect(or.crds, crs = crs)
por.crds <- project(or.crds, "epsg:4326")
o.cells <- cbind(o.cells, crds(por.crds))

ra.dat$cell <- o.cells[ , "cell"]
ra.dat$lon.centroid <- o.cells[ , "x"]
ra.dat$lat.centroid <- o.cells[ , "y"]

pal <- colorRampPalette(c("grey","white"))
plot(active.r, col = pal(2))
plot(proj.coords, add = TRUE)

# data source location map -----------------------------------------------------
library(sf)
library(rnaturalearth)

africa <- ne_countries(continent = 'Africa', scale = 'medium', returnclass = 'sf')

proj.coords.newcrs <- st_transform(st_as_sf(proj.coords), st_crs(africa))
same.crs(africa, proj.coords.newcrs)

library(ggplot2)
p <- ggplot(data = africa) +
  geom_sf() +
  theme_bw() +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  ggtitle("Location of relative abundance data sets") +
  geom_sf(
    data = proj.coords.newcrs,
    size = 1,
    shape = 16,
    colour = "blue"
  )
p
ggsave(
  "RelAbundDataMap.png",
  p,
  width = 20,
  height = 20,
  units = "cm"
)

# sample compilation -----------------------------------------------------------

# given aggregation to cells
# sample: same cell id, same date, same collection protocol, same origin year. Project agnostic.
ra.dat$samp.id <- paste(ra.dat$cell, ra.dat$Collection.protocols, ra.dat$origin.year)
samps <- sort(unique(ra.dat$samp.id))
# create dataframe of species frequencies
# roll through each sample
# preserve metadata from first record that appears in each sample
for (i in 1:length(samps)) {

  temp.ra <- ra.dat[ra.dat$samp.id == samps[i], ]
  temp.ra$Ag <- sum(temp.ra$Ag)
  temp.ra$Ac <- sum(temp.ra$Ac)
  temp.ra$Aa <- sum(temp.ra$Aa)

  # citation data and project id
  if(length(unique(temp.ra$Citations))>1) cat(i, " Non-unique citation \n")
  temp.ra$Citations <- paste(unique(temp.ra$Citations), collapse = " ; ")
  temp.ra$Projects <- paste(unique(temp.ra$Projects), collapse = " ; ")

  # null info
  temp.ra$Sample.ID <- temp.ra$Record.type <- temp.ra$Sample.type <-
    #temp.ra$Label <- temp.ra$Collection.date.range <- temp.ra$Projects <-
    temp.ra$Label <- temp.ra$Collection.date.range <-
    temp.ra$Latitudes <- temp.ra$Longitudes <- temp.ra$Locations <-
    #temp.ra$Citations <- temp.ra$Tag <- temp.ra$Attractants <-
    temp.ra$Tag <- temp.ra$Attractants <-
    temp.ra$Usage.license <- temp.ra$Sex <- temp.ra$Developmental.stage <-
    temp.ra$Available.data.types <- temp.ra$Geolocation.provenance <-
    temp.ra$Geolocation.precision <- NULL
  temp.ra <- temp.ra[1, ]

  if (i == 1) {
    # initialise abundance data set
    rel.abundance <- temp.ra[1, ]
  } else {
    rel.abundance <- rbind(rel.abundance, temp.ra)
  }

}

# Create a table of the citations for the relative abundance data --------------
## Many PMIDs associated with the same project
relAbund.cit <- unique(rel.abundance[, c(3, 4)])

## Identify the unique DOIs
relAbund.cit$DOI <- gsub(",.*$", "", relAbund.cit$Citations)
relAbund.cit <- unique(relAbund.cit[, c(1, 3)])

## Rename column headings
names(relAbund.cit) <- c("Project", "Citation")

## Export as a table
library(xtable)
print(
  xtable(relAbund.cit, caption = "Caption here"),
  include.rownames = F,
  size = "small",
  file = "RelAbundCit.txt"
)

# export relative abundance data -----------------------------------------------
rel.abundance <- rel.abundance[, -c(3,4)]
saveRDS(rel.abundance, file = "VB_PCR_rel_abundance_stra.rds")


