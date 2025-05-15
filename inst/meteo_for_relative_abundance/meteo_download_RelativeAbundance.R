# download of precipitations and temperature for relative abundance

# accessed 01/04/2025
# data are WGS84
# https://power.larc.nasa.gov/docs/methodology/

# https://power.larc.nasa.gov/#contact
# POWER data products are used in a publication, we request the follow-
#   ing acknowledgment be included: “These data were obtained from the NASA
# Langley Research Center POWER Project funded through the NASA Earth
# Science Directorate Applied Science Program.”

# relative abundance aggregated to grid cells
# use current ra file
ra <- readRDS("../01-RelativeAbundanceData/VB_PCR_rel_abundance_stra.rds")

# lat/long for each cell
meteo.cells <- unique(ra$cell)
meteo.lon <- meteo.lat <- mindate <- maxdate <- numeric(length(meteo.cells))

for (i in 1:length(meteo.cells)) {
  meteo.lon[i] <- unique(ra[ra$cell == meteo.cells[i], "lon.centroid"])
  meteo.lat[i] <- unique(ra[ra$cell == meteo.cells[i], "lat.centroid"])
  mindate[i] <- min(ra[ra$cell == meteo.cells[i], "date1"])
  maxdate[i] <- max(ra[ra$cell == meteo.cells[i], "date2"])
}
meteo.locations <- as.data.frame(cbind(meteo.cells, meteo.lon, meteo.lat,
                                        mindate,maxdate))

meteo.locations$meteo.lon <- as.numeric(meteo.locations$meteo.lon)
meteo.locations$meteo.lat <- as.numeric(meteo.locations$meteo.lat)

# within scope for NASA power
library(nasapower)

meteo.locations$PRECTOTCORR <- NA

pt <- proc.time()
for (i in 1:nrow(meteo.locations))
{
  # we are reading up to one year in the past for each day
  # in case we need to do a moving average over the last year
  start.date <- as.Date(meteo.locations$mindate[i]) - 365
  end.date <- meteo.locations$maxdate[i]

  clim.rec <- get_power(
    community = "ag",
    lonlat = c(meteo.locations$meteo.lon[i], meteo.locations$meteo.lat[i]),
    pars = c("T2M", "RH2M", "PRECTOTCORR"),
    dates = c(start.date, end.date),
    temporal_api = "daily"
  )

  # NASA/POWER CERES/MERRA2 Native Resolution Daily Data
  # The minimal and maximal dates depend of each location
  # Parameters:
  # T2M             MERRA-2 Temperature at 2 Meters (C) ;
  # RH2M            MERRA-2 Relative Humidity at 2 Meters (%) ;
  # PRECTOTCORR     MERRA-2 Precipitation Corrected (mm/day)

  write.csv(clim.rec,
            file = paste0("meteo_", meteo.cells[i], ".csv"))

  cat("completed ", i, " of ", nrow(meteo.locations), "\n")
  cat((proc.time()['elapsed'] - pt['elapsed'])/60, " min \n")

}


