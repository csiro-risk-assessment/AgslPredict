# spatial covariates -----------------------------------------------------------

library(terra)

# spatial covariates
# Lambert azimuthal equal area
# see Steinwand et al. (1995)

crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106
# read active grid
active.r <- rast("../covariates_spatial/activeAfrica.tif")
res(active.r) # 5000 5000 with skip = 1

# meteorological results -------------------------------------------------------

results <- readRDS("results.met.rds")
colnames(results) <- c(colnames(results[,1:14]),
                       c("PRECTOTCOR_Q1", "PRECTOTCOR_Q2", "PRECTOTCOR_Q3", "PRECTOTCOR_Q4"),
                       colnames(results[,19:33]))


# raster layers ----------------------------------------------------------------

results.r <- rast(active.r,
                  nlyrs = 29,
                  names = c(
                    "ActiveCell", "M2Cell",
                    paste("T2M_", "Q", 1:4, sep = ""),
                    paste("RH2M_", "Q", 1:4, sep = ""),
                    paste("PRECTOTCOR_", "Q", 1:4, sep = ""),
                    paste0(rep(c("T2M", "RH2M", "PRECTOTCOR"), each = 5),
                          c("_mean", paste0("_anom", seq(2002, 2020, 6))))
                  )
)

# cell id check
results.r[["ActiveCell"]][results[ , "cell"]] <- results[ , "cell"]
results.r[["ActiveCell"]][results[ , "cell"]] <- results[ , "M2CellNum"]

# quarterly means
results.r[["T2M_Q1"]][results[ , "cell"]] <- results[ , "T2M_Q1"]
results.r[["T2M_Q2"]][results[ , "cell"]] <- results[ , "T2M_Q2"]
results.r[["T2M_Q3"]][results[ , "cell"]] <- results[ , "T2M_Q3"]
results.r[["T2M_Q4"]][results[ , "cell"]] <- results[ , "T2M_Q4"]
results.r[["RH2M_Q1"]][results[ , "cell"]] <- results[ , "RH2M_Q1"]
results.r[["RH2M_Q2"]][results[ , "cell"]] <- results[ , "RH2M_Q2"]
results.r[["RH2M_Q3"]][results[ , "cell"]] <- results[ , "RH2M_Q3"]
results.r[["RH2M_Q4"]][results[ , "cell"]] <- results[ , "RH2M_Q4"]
results.r[["PRECTOTCOR_Q1"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_Q1"]
results.r[["PRECTOTCOR_Q2"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_Q2"]
results.r[["PRECTOTCOR_Q3"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_Q3"]
results.r[["PRECTOTCOR_Q4"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_Q4"]
# annual means
results.r[["T2M_mean"]][results[ , "cell"]] <- results[ , "T2M_mean"]
results.r[["RH2M_mean"]][results[ , "cell"]] <- results[ , "RH2M_mean"]
results.r[["PRECTOTCOR_mean"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_mean"]
# anomalies
results.r[["T2M_anom2002"]][results[ , "cell"]] <- results[ , "T2M_y2002"] - results[ , "T2M_mean"]
results.r[["T2M_anom2008"]][results[ , "cell"]] <- results[ , "T2M_y2008"] - results[ , "T2M_mean"]
results.r[["T2M_anom2014"]][results[ , "cell"]] <- results[ , "T2M_y2014"] - results[ , "T2M_mean"]
results.r[["T2M_anom2020"]][results[ , "cell"]] <- results[ , "T2M_y2020"] - results[ , "T2M_mean"]
results.r[["RH2M_anom2002"]][results[ , "cell"]] <- results[ , "RH2M_y2002"] - results[ , "RH2M_mean"]
results.r[["RH2M_anom2008"]][results[ , "cell"]] <- results[ , "RH2M_y2008"] - results[ , "RH2M_mean"]
results.r[["RH2M_anom2014"]][results[ , "cell"]] <- results[ , "RH2M_y2014"] - results[ , "RH2M_mean"]
results.r[["RH2M_anom2020"]][results[ , "cell"]] <- results[ , "RH2M_y2020"] - results[ , "RH2M_mean"]
results.r[["PRECTOTCOR_anom2002"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_y2002"] - results[ , "PRECTOTCOR_mean"]
results.r[["PRECTOTCOR_anom2008"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_y2008"] - results[ , "PRECTOTCOR_mean"]
results.r[["PRECTOTCOR_anom2014"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_y2014"] - results[ , "PRECTOTCOR_mean"]
results.r[["PRECTOTCOR_anom2020"]][results[ , "cell"]] <- results[ , "PRECTOTCOR_y2020"] - results[ , "PRECTOTCOR_mean"]

# projection to WGS84
res.r <- project(results.r, "epsg:4326")

# plots ------------------------------------------------------------------------

library(rasterVis)
#
# RColorBrewer::display.brewer.all(n=7, type="seq", select=NULL, exact.n=TRUE,
#                                                      colorblindFriendly=TRUE)
# RColorBrewer::display.brewer.all(n=7, type="div", select=NULL, exact.n=TRUE,
#                                  colorblindFriendly=TRUE)

## annual ----------------------------------------------------------------------

brks.RH2M <- c(10*(0:10))
tcks.RH2M <- 0:10
labs.RH2M <- c(expression("0%"), expression("10%"), expression("20%"),
               expression("30%"), expression("40%"), expression("50%"),
               expression("60%"), expression("70%"), expression("80%"),
               expression("90%"), expression("100%")
          )

brks.T2M <- c(5*(0:7))
tcks.T2M <- 0:7
labs.T2M <- c(expression("0"), expression("5"),
              expression("10"), expression("15"),
              expression("20"), expression("25"),
              expression("30"), expression("35")
              )

brks.PR <- c(0:10)
tcks.PR <- 0:10
labs.PR <- c(expression("0"), expression("1"), expression("2"),
             expression("3"), expression("4"), expression("5"),
             expression("6"), expression("7"), expression("8"),
             expression("9"), expression("10")
)

cols <- c(RColorBrewer::brewer.pal(9, "PuBu"))
rtheme <- rasterTheme(region = cols)
ann.theme.RH2M <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

cols <- c(RColorBrewer::brewer.pal(9, "YlOrRd"))
rtheme <- rasterTheme(region = cols)
ann.theme.T2M <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

cols <- c(RColorBrewer::brewer.pal(9, "Blues"))
rtheme <- rasterTheme(region = cols)
ann.theme.PR <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)

png("annual_avg_T2M.png", 2.5*480, 2.5*480*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("T2M_mean")]],
          main =  "Annual average temperature",
          par.settings = ann.theme.T2M,
          at = brks.T2M,
          colorkey = list(
            at = tcks.T2M,
            labels = list(
              at = tcks.T2M,
              labels = labs.T2M
            )
          ),
          margin = FALSE
)
dev.off()

png("annual_avg_RH2M.png", 2.5*480, 2.5*480*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("RH2M_mean")]],
          main = "Annual average relative humidity",
          par.settings = ann.theme.RH2M,
          at = brks.RH2M,
          colorkey = list(
            at = tcks.RH2M,
            labels = list(
              at = tcks.RH2M,
              labels = labs.RH2M
            )
          ),
          margin = FALSE
)
dev.off()

png("annual_avg_PRECTOTCOR.png", 2.5*480, 2.5*480*(0.8), pointsize = 72, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("PRECTOTCOR_mean")]],
     main = "Annual average precipitation",
     par.settings = ann.theme.PR,
     at = brks.PR,
     colorkey = list(
       at = tcks.PR,
       labels = list(
         at = tcks.PR,
         labels = labs.PR
       )
     ),
     margin = FALSE
)
dev.off()

## quarterly -------------------------------------------------------------------
cols <- c(RColorBrewer::brewer.pal(11, "PuBu"))
quart.theme.RH2M <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)

cols <- c(RColorBrewer::brewer.pal(11, "YlOrRd"))
quart.theme.T2M <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)

cols <- c(RColorBrewer::brewer.pal(11, "Blues"))
quart.theme.PR <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)

brks.PRQ <- seq(from = 0, to = 20, length.out = 11)
tcks.PRQ <- seq(from = 0, to = 20, length.out = 11)
labs.PRQ <- c(expression("0"), expression("2"), expression("4"),
             expression("6"), expression("8"), expression("10"),
             expression("12"), expression("14"), expression("16"),
             expression("18"), expression("20")
)

png("Q_T2M.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("T2M_Q1", "T2M_Q2", "T2M_Q3", "T2M_Q4")]],
          main = expression("Quarterly average temperature"),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme.T2M,
          at = brks.T2M,
          colorkey = list(
            at = tcks.T2M,
            labels = list(
              at = tcks.T2M,
              labels = labs.T2M
            )
          )
)
dev.off()
png("Q_RH2M.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("RH2M_Q1", "RH2M_Q2", "RH2M_Q3", "RH2M_Q4")]],
          main = expression("Quarterly average relative humidity"),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme.RH2M,
          at = brks.RH2M,
          colorkey = list(
            at = tcks.RH2M,
            labels = list(
              at = tcks.RH2M,
              labels = labs.RH2M
            )
          )
)
dev.off()
png("Q_PRECTOTCOR.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("PRECTOTCOR_Q1", "PRECTOTCOR_Q2", "PRECTOTCOR_Q3", "PRECTOTCOR_Q4")]],
          main = expression("Quarterly average precipitation"),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme.PR,
          at = brks.PRQ,
          colorkey = list(
            at = tcks.PRQ,
            labels = list(
              at = tcks.PRQ,
              labels = labs.PRQ
            )
          )
)
dev.off()

## anomalies -------------------------------------------------------------------

brks.RH2M <- c(seq(-20,20,4))
tcks.RH2M <- c(seq(-20,20,4))
labs.RH2M <- c(expression("-20%"), expression("-16%"), expression("-12%"),
               expression("-8%"), expression("-4%"), expression("0%"),
               expression("4%"), expression("8%"), expression("12%"),
               expression("16%"), expression("20%")
)

brks.T2M <- c(5e-1*(-5:5))
tcks.T2M <- c(5e-1*(-5:5))
labs.T2M <- c(expression("-2.5"),expression("-2.0"),
              expression("-1.5"),expression("-1.0"),
              expression("-0.5"),expression("0.0"),
              expression("0.5"), expression("1.0"),
              expression("1.5"), expression("2.0"),
              expression("2.5")
)

brks.PR <- c(seq(-15,15,3))
tcks.PR <- c(seq(-15,15,3))
labs.PR <- c(expression("-15"), expression("-12"), expression("-9"),
             expression("-6"), expression("-3"), expression("0"),
             expression("3"), expression("6"), expression("9"),
             expression("12"), expression("15")
)

cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
anom.theme <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)
png("Anom_T2M.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("T2M_anom2002", "T2M_anom2008", "T2M_anom2014", "T2M_anom2020")]],
          main = expression("Anomaly relative to temperature 2002-2020 average"),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks.T2M,
          colorkey = list(
            at = tcks.T2M,
            labels = list(
              at = tcks.T2M,
              labels = labs.T2M
            )
          )
)
dev.off()
png("Anom_RH2M.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("RH2M_anom2002", "RH2M_anom2008", "RH2M_anom2014", "RH2M_anom2020")]],
          main = expression("Anomaly relative to 2002-2020 relative humidity average"),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks.RH2M,
          colorkey = list(
            at = tcks.RH2M,
            labels = list(
              at = tcks.RH2M,
              labels = labs.RH2M
            )
          )
)
dev.off()
png("Anom_PRECTOTCOR.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("PRECTOTCOR_anom2002", "PRECTOTCOR_anom2008", "PRECTOTCOR_anom2014", "PRECTOTCOR_anom2020")]],
          main = expression("Anomaly relative to 2002-2020 precipitation average"),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks.PR,
          colorkey = list(
            at = tcks.PR,
            labels = list(
              at = tcks.PR,
              labels = labs.PR
            )
          )
)
dev.off()

