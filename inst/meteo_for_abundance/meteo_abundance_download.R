# download of precip for abundance analysis

# accessed 2 April 2025
# data are WGS84
# https://power.larc.nasa.gov/docs/methodology/

# https://power.larc.nasa.gov/#contact
# POWER data products are used in a publication, we request the follow-
#   ing acknowledgment be included: “These data were obtained from the NASA
# Langley Research Center POWER Project funded through the NASA Earth
# Science Directorate Applied Science Program.”

# abundance aggregated to grid cells
abund <- readRDS("../03-AbundanceData/abundance.rds")

# duration "1997-01-14" to "2021-04-27"
# within scope for NASA power
start.date <- "1996-07-01"
end.date <- "2021-04-27"

# lat/long for each cell
precip.cells <- unique(abund$cell)
precip.lon <- precip.lat <- numeric(length(precip.cells))
for (i in 1:length(precip.cells)) {
  precip.lon[i] <- unique(abund[abund$cell == precip.cells[i], "lon.centroid"])
  precip.lat[i] <- unique(abund[abund$cell == precip.cells[i], "lat.centroid"])
}
precip.locations <- as.data.frame(cbind(precip.cells, precip.lon, precip.lat))

library(nasapower) # v4.2.1

query_parameters(community = "ag", par = "PRECTOTCORR", temporal_api ="daily")

precip.locations$PRECTOTCORR <- NA

pt <- proc.time()
for (i in 1:nrow(precip.locations)) {
  
  clim.rec <- get_power(
    community = "ag",
    lonlat = c(precip.locations$precip.lon[i], precip.locations$precip.lat[i]),
    pars = c("T2M", "RH2M", "PRECTOTCORR"),
    dates = c(start.date, end.date),
    temporal_api = "daily"
  )
  
  # NASA/POWER CERES/MERRA2 Native Resolution Daily Data  
  # Dates (month/day/year): 07/01/1996 through 04/27/2021 
  # Parameters: 
  # T2M             MERRA-2 Temperature at 2 Meters (C) ;
  # RH2M            MERRA-2 Relative Humidity at 2 Meters (%) ;
  #   PRECTOTCORR     MERRA-2 Precipitation Corrected (mm/day)
  
  write.csv(clim.rec, 
            file = paste0("precip_", precip.cells[i], ".csv"))
  
  cat("completed ", i, " of ", nrow(precip.locations), "\n")
  cat((proc.time()['elapsed'] - pt['elapsed'])/60, " min \n")
  
}


