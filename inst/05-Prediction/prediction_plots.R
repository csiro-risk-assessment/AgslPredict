# prediction plots
# using results from prediction_tot_relabund_array.R
# and tempMask.R

# read in predictions ----------------------------------------------------------

# source array of results and conatenate into single results file
res.files <- grep("results_", dir(), value = TRUE)
results <- readRDS(res.files[1])
results[ , 5:ncol(results)] <- NA
for (i in 1:length(res.files)) {
  res <- readRDS(res.files[i])
  res.ids <- which(!is.na(res[ , 6]))
  if (any(!is.na(results[res.ids, 5:ncol(results)])))
    stop("indexing error")
  if (!all(results[res.ids, 1:4] == res[res.ids, 1:4]))
    stop("results structure error")
  results[res.ids, 5:ncol(res)] <- res[res.ids, 5:ncol(res)] 
}

# human population -------------------------------------------------------------
crs <- "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84" # EPSG:42106

library(terra)

active.r <- rast("../covariates_spatial/activeAfrica.tif")
same.crs(crs(active.r), crs) # TRUE
res(active.r) 

## read in remaining spatial covariate .csv files

# check for correct extent info
all(read.csv("../covariates_spatial/pop.csv", header = FALSE, nrow = 1) ==
      c("#xmin=-4099134.0",	"ymin=-4202349.0", "cell_size=5000.0",	"nx=1520", "ny=1280"))
pop <- read.csv("../covariates_spatial/pop.csv", skip = 1, header = FALSE)
all.equal(dim(pop), c(1280, 1520))
pop <- as.matrix(pop)
pop.r <- rast(pop[nrow(pop):1, ], crs = crs,
              extent = c(-4099134.0, -4099134.0+5000*1520, -4202349.0, -4202349.0+5000*1280))

res(pop.r)

# temperature mask -------------------------------------------------------------

tmask <- readRDS("threshmat.rds")
results[tmask[ , "TOK"] == 0, c(
  "Aa_Q1",     "Aa_Q2",     "Aa_Q3",     "Aa_Q4",     "Ac_Q1",     "Ac_Q2",
  "Ac_Q3",     "Ac_Q4",     "Ag_Q1",     "Ag_Q2",     "Ag_Q3",     "Ag_Q4",
  "Aa_mean",   "Aa_y2002",  "Aa_y2008",  "Aa_y2014",  "Aa_y2020",  "Ac_mean" ,
  "Ac_y2002",  "Ac_y2008",  "Ac_y2014",  "Ac_y2020",  "Ag_mean",   "Ag_y2002",
  "Ag_y2008",  "Ag_y2014",  "Ag_y2020"
)] <- 0

# total pop --------------------------------------------------------------------

# A recommendation on the method to delineate cities, urban and rural areas
# for international statistical comparisons
# Prepared by the European Commission – Eurostat and DG for Regional and Urban Policy –
# ILO, FAO, OECD, UN-Habitat, World Bank
urban.threshold <- 1500 # per km^2
sum(values(pop.r) > urban.threshold, na.rm = TRUE)/sum(!is.na(values(pop.r)))

thresh.ids <- which(values(pop.r) > 1500)

# mask urban areas
results.nu <- results
results.nu[results.nu[ , "cell"]%in%thresh.ids, 7:ncol(results.nu)] <- NA

formatC(sum(results.nu[ , "Aa_mean"], na.rm = TRUE))
formatC(sum(results.nu[ , "Ac_mean"], na.rm = TRUE))
formatC(sum(results.nu[ , "Ag_mean"], na.rm = TRUE))

# adult population size estimates: both sexes
# Aa, Ac, Ag by year and overall average
# Ac + Ag effective population size
pop.tab <- matrix(NA, nrow = 5, ncol = 1,
                  dimnames = list(Years = c("2002", "2008", "2014", "2020", "Mean"),
                                  Population = c("Effective")))
# Ac + Ag
# adjust for rough ratio of effective population size to total Ac + Ag (female + male)
pop.tab["2002", c("Effective")] <-
  formatC(sum(rowSums(results.nu[ , c("Ac_y2002", "Ag_y2002")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2008", c("Effective")] <-
  formatC(sum(rowSums(results.nu[ , c("Ac_y2008", "Ag_y2008")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2014", c("Effective")] <-
  formatC(sum(rowSums(results.nu[ , c("Ac_y2014", "Ag_y2014")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2020", c("Effective")] <-
  formatC(sum(rowSums(results.nu[ , c("Ac_y2020", "Ag_y2020")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["Mean", c("Effective")] <-
  formatC(sum(rowSums(results.nu[ , c("Ac_mean", "Ag_mean")]*2*0.10, na.rm = TRUE)), digits = 3)

library(xtable)

s <- sanitize.numbers(pop.tab, type = "latex", math.style.exponents = TRUE)

# EIR
hpop <- extract(pop.r, results[ , c('x', 'y')])

# approx species feeding rates from daily feeding probabilities
rate <- matrix(-log((1 - c(0.37, 0.49, 0.47))), nrow = nrow(results.nu), ncol = 3, byrow = TRUE)

# sporozoite rates for An. gambiae s.l. from Table 3 in Msugupakulya et al. (2023)
bv2002 <- rowSums(results.nu[ , c("Aa_y2002", "Ac_y2002", "Ag_y2002")]*rate, na.rm = TRUE)*0.015 
bv2020 <- rowSums(results.nu[ , c("Aa_y2020", "Ac_y2020", "Ag_y2020")]*rate, na.rm = TRUE)*0.010
# divide by human population
bv2002 <- bv2002/(hpop$lyr.1*25)
bv2020 <- bv2020/(hpop$lyr.1*25)
# annualise 
bv2002 <- bv2002*365
bv2020 <- bv2020*365
# # # filter to at least 50 inhabitants per sq km
# bv2002 <- bv2002[hp]
# bv2020 <- bv2020[hp]

qs <- c(0.1, 0.5, 0.9)
quantile(bv2002[results.nu[ , "lon"] >= 30], qs, na.rm = TRUE)
quantile(bv2020[results.nu[ , "lon"] >= 30], qs, na.rm = TRUE)

saveRDS(results.nu, file = "resultsNU.rds")

# raster layers ----------------------------------------------------------------

results.r <- rast(active.r,
                  nlyrs = 36,
                  names = c(
                    "activeCell", "M2Cell",
                    paste("Aa_", "Q", 1:4, sep = ""),
                    paste("Ac_", "Q", 1:4, sep = ""),
                    paste("Ag_", "Q", 1:4, sep = ""),
                    paste0(rep(c("Aa", "Ac", "Ag"), each = 5),
                          c("_mean", paste0("_anom", seq(2002, 2020, 6)))),
                    paste0(c("Aa", "Ac", "Ag"), "_r"),
                    c("logAcAa_v", "logAgAa_v", "tot_mean", "tot_v")
                  )
)

# cell id check
results.r[["activeCell"]][results.nu[ , "cell"]] <- results.nu[ , "cell"]
results.r[["activeCell"]][results.nu[ , "cell"]] <- results.nu[ , "M2CellNum"]

# quarterly means
results.r[["Aa_Q1"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_Q1"]
results.r[["Aa_Q2"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_Q2"]
results.r[["Aa_Q3"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_Q3"]
results.r[["Aa_Q4"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_Q4"]
results.r[["Ac_Q1"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_Q1"]
results.r[["Ac_Q2"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_Q2"]
results.r[["Ac_Q3"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_Q3"]
results.r[["Ac_Q4"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_Q4"]
results.r[["Ag_Q1"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_Q1"]
results.r[["Ag_Q2"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_Q2"]
results.r[["Ag_Q3"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_Q3"]
results.r[["Ag_Q4"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_Q4"]
# annual means
results.r[["Aa_mean"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_mean"]
results.r[["Ac_mean"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_mean"]
results.r[["Ag_mean"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_mean"]
# anomalies
results.r[["Aa_anom2002"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_y2002"] - results.nu[ , "Aa_mean"]
results.r[["Aa_anom2008"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_y2008"] - results.nu[ , "Aa_mean"]
results.r[["Aa_anom2014"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_y2014"] - results.nu[ , "Aa_mean"]
results.r[["Aa_anom2020"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_y2020"] - results.nu[ , "Aa_mean"]
results.r[["Ac_anom2002"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_y2002"] - results.nu[ , "Ac_mean"]
results.r[["Ac_anom2008"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_y2008"] - results.nu[ , "Ac_mean"]
results.r[["Ac_anom2014"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_y2014"] - results.nu[ , "Ac_mean"]
results.r[["Ac_anom2020"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_y2020"] - results.nu[ , "Ac_mean"]
results.r[["Ag_anom2002"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_y2002"] - results.nu[ , "Ag_mean"]
results.r[["Ag_anom2008"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_y2008"] - results.nu[ , "Ag_mean"]
results.r[["Ag_anom2014"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_y2014"] - results.nu[ , "Ag_mean"]
results.r[["Ag_anom2020"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_y2020"] - results.nu[ , "Ag_mean"]
# raw predictions
results.r[["Aa_r"]][results.nu[ , "cell"]] <- results.nu[ , "Aa_r"]
results.r[["Ac_r"]][results.nu[ , "cell"]] <- results.nu[ , "Ac_r"]
results.r[["Ag_r"]][results.nu[ , "cell"]] <- results.nu[ , "Ag_r"]
results.r[["logAcAa_v"]][results.nu[ , "cell"]] <- results.nu[ , "logAcAa_v"]
results.r[["logAgAa_v"]][results.nu[ , "cell"]] <- results.nu[ , "logAgAa_v"]
results.r[["tot_mean"]][results.nu[ , "cell"]] <- results.nu[ , "tot_mean"]
results.r[["tot_v"]][results.nu[ , "cell"]] <- results.nu[ , "tot_v"]

# projection to WGS84
res.r <- project(results.r, "epsg:4326")

# plots ------------------------------------------------------------------------

library(rnaturalearth)
library(rasterVis)
library(latticeExtra)

RColorBrewer::display.brewer.all(n=7, type="seq", select=NULL, exact.n=TRUE,
                                                     colorblindFriendly=TRUE)
RColorBrewer::display.brewer.all(n=7, type="div", select=NULL, exact.n=TRUE,
                                 colorblindFriendly=TRUE)

# shoreline from rnaturalearth
shore <- ne_coastline(scale = "medium", returnclass = "sv")

## annual ----------------------------------------------------------------------

brks <- c(0, 10^(0:7))
tcks <- 0:8
labs <- c(0, 1, 10,
          expression(10^2), expression(10^3), expression(10^4),
          expression(10^5), expression(10^6), expression(10^7)
          )
cols <- c(grey(0.5), RColorBrewer::brewer.pal(7, "YlOrRd"))
rtheme <- rasterTheme(region = cols)
ann.theme <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)
png("annual_avg.png", 2.5*480*3, 2.5*480*(1/3)*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("Aa_mean", "Ac_mean", "Ag_mean")]],
          maxpixels = 2.3e6,
          main = "Annual average abundance",
          names.attr = c(
            expression(italic("An. arabiensis")),
            expression(italic("An. coluzzii")),
            expression(italic("An. gambiae")~s.s.)
          ),
          layout = c(3, 1),
          par.settings = ann.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()

## quarterly -------------------------------------------------------------------

quart.theme <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)
png("Q_Aa.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Aa_Q1", "Aa_Q2", "Aa_Q3", "Aa_Q4")]],
          maxpixel = 2.3e6,
          main = expression(paste("Quarterly average abundance: ",
                                  italic("An. arabiensis"))),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()
png("Q_Ac.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Ac_Q1", "Ac_Q2", "Ac_Q3", "Ac_Q4")]],
          maxpixel = 2.3e6,
          main = expression(paste("Quarterly average abundance: ",
                                  italic("An. coluzzii"))),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()
png("Q_Ag.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Ag_Q1", "Ag_Q2", "Ag_Q3", "Ag_Q4")]],
          maxpixel = 2.3e6,
          main = expression(paste("Quarterly average abundance: ",
                                  italic("An. gambiae"), " s.s.")),
          names.attr = paste0("Q", 1:4),
          layout = c(2, 2),
          par.settings = quart.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()

## anomalies -------------------------------------------------------------------

brks <- c(-10^(seq(8, 2, -2)), 0, 10^(seq(2, 8, 2)))
tcks <- c(-seq(8, 2, -2), 0, seq(2, 8, 2))
labs <- c(
  expression(-10^8), expression(-10^6),
  expression(-10^4), expression(-10^2),
  0,
  expression(10^2), expression(10^4),
  expression(10^6), expression(10^8)
)
cols <- rev(RColorBrewer::brewer.pal(8, "RdBu"))
anom.theme <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2),
  fontsize = list(
    text = 18, points = 8
  )
)
png("Anom_Aa.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Aa_anom2002", "Aa_anom2008", "Aa_anom2014", "Aa_anom2020")]],
          maxpixel = 2.3e6,
          main = expression(paste("Anomaly relative to 2002-2020 average: ",
                                  italic("An. arabiensis"))),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()
png("Anom_Ac.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Ac_anom2002", "Ac_anom2008", "Ac_anom2014", "Ac_anom2020")]],
          maxpixel = 2.3e6,
          main = expression(paste("Anomaly relative to 2002-2020 average: ",
                                  italic("An. coluzzii"))),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()
png("Anom_Ag.png", 2.5*480*3, 2.5*480*(0.8)*3, pointsize = 72, res = 3*72)
levelplot(res.r[[c("Ag_anom2002", "Ag_anom2008", "Ag_anom2014", "Ag_anom2020")]],
          maxpixel = 2.3e6,
          main = expression(paste("Anomaly relative to 2002-2020 average: ",
                                  italic("An. gambiae"), " s.s.")),
          names.attr = c("2002", "2008", "2014", "2020"),
          layout = c(2, 2),
          par.settings = anom.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()

## raw rel abundance------------------------------------------------------------

brks <- seq(from = 0, to = 1, length.out = 9)
tcks <- brks
labs <- brks
cols <- c(grey(0.5), RColorBrewer::brewer.pal(7, "YlOrRd"))
rtheme <- rasterTheme(region = cols)
ann.theme <- rasterTheme(
  region = cols,
  layout.heights=list(top.padding = -2,
                      bottom.padding = -2)
)
png("rel_abundance.png", 2.5*480*3, 2.5*480*(1/3)*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("Aa_r", "Ac_r", "Ag_r")]],
          maxpixel = 2.3e6,
          main = "Predicted relative abundance",
          names.attr = c(
            expression(italic("An. arabiensis")),
            expression(italic("An. coluzzii")),
            expression(italic("An. gambiae")~s.s.)
          ),
          layout = c(3, 1),
          par.settings = ann.theme,
          at = brks,
          colorkey = list(
            at = tcks,
            labels = list(
              at = tcks,
              labels = labs
            )
          )
) + layer(llines(shore, col = "black"))
dev.off()

## predicted variance raw rel abundance ----------------------------------------

png("rel_abundance_var.png", 2*480*3, 2*480*(1/2)*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[[c("logAcAa_v", "logAgAa_v")]],
          maxpixel = 2.3e6,
          main = "Average prediction variance",
          names.attr = c(
            expression(paste("log ", italic("An. coluzzii") / italic("An. arabiensis"))),
            expression(paste("log ", italic("An. gambiae")~s.s. / italic("An. arabiensis")))
          ),
          layout = c(2, 1),
          zscaleLog = TRUE
) + layer(llines(shore, col = "black"))
dev.off()


## raw total abundance ---------------------------------------------------------

png("tot_abundance.png", sqrt(0.5)*2*480*3, sqrt(0.5)*2*480*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(log(res.r[["tot_mean"]]),
          maxpixel = 2.3e6,
          main = "Predicted per person average abundance",
          margin = FALSE,
          zscaleLog = FALSE,
          colorkey = list(at = log(seq(1, 32, 0.5)),
                        labels = list(
                          at = log(c(1, 2, 4, 8, 16, 32)),
                          labels = c(1, 2, 4, 8, 16, 32)
                          )
          )
) + layer(llines(shore, col = "black"))
dev.off()

## predicted variance total abundance ------------------------------------------

png("tot_abundance_var.png", sqrt(0.5)*2*480*3, sqrt(0.5)*2*480*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[["tot_v"]],
          maxpixel = 2.3e6,
          main = "Predicted per person average variance of log abundance",
          margin = FALSE,
          zscaleLog = TRUE
) + layer(llines(shore, col = "black"))
dev.off()

