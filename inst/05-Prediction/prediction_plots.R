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

# active cells -----------------------------------------------------------------

library(terra)
active.r <- rast("../covariates_spatial/activeAfrica.tif")

# temperature mask -------------------------------------------------------------
# produced by script "tempMask.R"

tmask <- readRDS("threshmat.rds")
results[tmask[ , "TOK"] == 0, c(
  "Aa_Q1",     "Aa_Q2",     "Aa_Q3",     "Aa_Q4",     "Ac_Q1",     "Ac_Q2",
  "Ac_Q3",     "Ac_Q4",     "Ag_Q1",     "Ag_Q2",     "Ag_Q3",     "Ag_Q4",
  "Aa_mean",   "Aa_y2002",  "Aa_y2008",  "Aa_y2014",  "Aa_y2020",  "Ac_mean" ,
  "Ac_y2002",  "Ac_y2008",  "Ac_y2014",  "Ac_y2020",  "Ag_mean",   "Ag_y2002",
  "Ag_y2008",  "Ag_y2014",  "Ag_y2020"
)] <- 0

# adult effective population size estimates: both sexes -------------------------
# Ac + Ag effective population size
pop.tab <- matrix(NA, nrow = 4, ncol = 1,
                  dimnames = list(Years = c("2002", "2008", "2014", "2020"),
                                  Population = c("Effective")))
# Ac + Ag
# adjust for rough ratio of effective population size to total Ac + Ag (female + male)
# thresholding high density urban cells
pop.tab["2002", c("Effective")] <-
  formatC(sum(rowSums(results[ , c("Ac_y2002", "Ag_y2002")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2008", c("Effective")] <-
  formatC(sum(rowSums(results[ , c("Ac_y2008", "Ag_y2008")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2014", c("Effective")] <-
  formatC(sum(rowSums(results[ , c("Ac_y2014", "Ag_y2014")]*2*0.10, na.rm = TRUE)), digits = 3)
pop.tab["2020", c("Effective")] <-
  formatC(sum(rowSums(results[ , c("Ac_y2020", "Ag_y2020")]*2*0.10, na.rm = TRUE)), digits = 3)

library(xtable)

s <- sanitize.numbers(pop.tab, type = "latex", math.style.exponents = TRUE)
print(xtable(s), file = "AcAg_effective.txt")

# EIR comparison East Africa ---------------------------------------------------
pop.files <- paste0("../human_pop/human_pop_", 2002:2020, ".tif")
pop.r <- rast(pop.files)
names(pop.r) <- paste0("Y", 2002:2020)
hpop <- extract(pop.r, results[ , c('x', 'y')])

# approx species feeding rates from daily feeding probabilities
rate <- matrix(-log((1 - c(0.37, 0.49, 0.47))), nrow = nrow(results), ncol = 3, byrow = TRUE)

# sporozoite rates for An. gambiae s.l. from Table 3 in Msugupakulya et al. (2023)
bv2002 <- rowSums(results[ , c("Aa_y2002", "Ac_y2002", "Ag_y2002")]*rate, na.rm = TRUE)*0.015 
bv2020 <- rowSums(results[ , c("Aa_y2020", "Ac_y2020", "Ag_y2020")]*rate, na.rm = TRUE)*0.010
# divide by human population --> cells with zero population NaN for EIR comparison
bv2002 <- bv2002/(hpop["Y2002"])
bv2020 <- bv2020/(hpop["Y2020"])
# annualise 
bv2002 <- bv2002*365
bv2020 <- bv2020*365

qs <- c(0.1, 0.5, 0.9)
EIR.tab <- matrix(NA, nrow = 3, ncol = 2,
                  dimnames = list(q = qs,
                                  EIR = c("2002",
                                    "2020" 
                                  ))
                  )
# window: east of 28 degrees E and south of 4 degrees N
# Fig. 5 in Msugupakulya et al. (2023)
EIR.tab[ , "2002"] <-
  quantile(bv2002[results[ , "lon"] >= 28 & results[ , "lat"] <= 4, ], qs, na.rm = TRUE)
EIR.tab[ , "2020"] <-
  quantile(bv2020[results[ , "lon"] >= 28 & results[ , "lat"] <= 4, ], qs, na.rm = TRUE)
print(xtable(EIR.tab, digits = 1), file = "EIR.txt")

# results for spatial plots
results.nu <- results
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

# shoreline from rnaturalearth
shore <- ne_coastline(scale = "medium", returnclass = "sv")

## annual ----------------------------------------------------------------------

brks <- c(0, 10^(0:8))
tcks <- 0:9
labs <- c(0, 1, 10,
          expression(10^2), expression(10^3), expression(10^4),
          expression(10^5), expression(10^6), expression(10^7),
          expression(10^8)
          )
cols <- c(grey(0.5), RColorBrewer::brewer.pal(8, "YlOrRd"))
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
          main = "Average daily abundance",
          names.attr = c(
            expression(italic("An. arabiensis")),
            expression(italic("An. coluzzii")),
            expression(italic("An. gambiae s.s."))
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
          main = expression(paste("Average daily abundance by quarter: ",
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
          main = expression(paste("Average daily abundance by quarter: ",
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
          main = expression(paste("Average daily abundance by quarter: ",
                                  italic("An. gambiae s.s."))),
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
                                  italic("An. gambiae s.s."))),
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

ann.theme <- rasterTheme(
  region = viridis::viridis(20),
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
            expression(italic("An. gambiae s.s."))
          ),
          layout = c(3, 1),
          par.settings = ann.theme,
          # at = brks,
          colorkey = list(at = seq(0, 1, 0.05),
                          labels = list(
                            at = seq(0, 1, by = 0.1),
                            labels = seq(0, 1, by = 0.1)
                          ))
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
            expression(paste("log ", italic("An. gambiae s.s.") / italic("An. arabiensis")))
          ),
          layout = c(2, 1),
          zscaleLog = TRUE
) + layer(llines(shore, col = "black"))
dev.off()


## raw total abundance ---------------------------------------------------------

png("tot_abundance.png", sqrt(0.5)*2*480*3, sqrt(0.5)*2*480*3, res = 3*72)
par(mar = rep(0.1, 4), oma = rep(1, 4))
levelplot(res.r[["tot_mean"]],
          maxpixel = 2.3e6,
          main = "Predicted per person average abundance",
          margin = FALSE,
          zscaleLog = TRUE
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

