# relative abundance analysis
# multinomial model
# one month only analysis

dat1m <- readRDS("../01-RelativeAbundanceData/ra1m_obs_covs_stra.rds")

# Begin with comparison to previous analysis: a time invariant covariate example.
# In this section these data are analysed with the model structure proposed for
# the multinomial GLM analysis of relative abundance in Appendix C.5 of:
# Hosack, G. R., Beeton, N. J., Ickowicz, A., Peel, D., Wilkins, A., Dambacher, J. M.,
# Wickramarachchi, A., McDonald, M., Tay, W. T., Wilson, L., Bauer, D., and Hayes, K.
# R. 2023. Risk Assessment for Controlling Mosquito Vectors with Engineered Nucleases:
#   Paternal Male Bias Construct. Report No. EP2022-4945. Hobart, Australia: CSIRO.
# https://doi.org/10.25919/2t8h-5k81
# Note: For the proposed model in Hosack et al. (2023) the VectorBase data were
# downloaded in August 2023 using a different user interface and structure as
# was available on VectorBase.org at that time. These data were not filtered by
# duration of "sample". The comparative analysis is included here for comparison
# purposes only and  applies to the data downloaded to directory
# system.file(package = "AgslPredict", "01-RelativeAbundanceData/")
# that occurred on 24 April 2024 using the archived version of the legacy map
# view dated to July 2023 then available on VectorBase.

# time invariant (historical) comparison ---------------------------------------

library(nnet)

# for this "historical" analysis use max 3 month duration of sample
dat3m <- readRDS("../01-RelativeAbundanceData/ra3m_obs_covs_stra.rds")
histDat <- dat3m

# frequency counts of species
histCounts <- cbind(histDat$Aa, histDat$Ac, histDat$Ag)

# predictors
histPreds <- histDat[ , c("lon.centroid", "lat.centroid",
                          "d2c", "d2r", "pop", "elev", "origin.year")]

# TRANSFORMATIONS
# log1p transform for all spatial covariates (except geographic coords)
histPreds[ , c("d2c", "d2r", "pop", "elev")] <-
  log1p(histPreds[ , c("d2c", "d2r", "pop", "elev")])
# absolute value of latitude (relating to insolation)
histPreds[ , "lat.centroid"] <- abs(histPreds[ , "lat.centroid"])

# fit model
histMod <- multinom(histCounts ~ (lon.centroid + lat.centroid + d2r + pop + elev)^2 +
                    + I(lon.centroid^2) + I(lat.centroid^2) + I(d2r^2) + I(pop^2) + I(elev^2) +
                    + origin.year, data = histPreds, maxit = 200)

# save model estimates
histEst <- cbind(t(summary(histMod)$coefficients), t(summary(histMod)$standard.errors))
colnames(histEst) <- c("coeffAcVsAa", "coeffAgVsAa", "seAcVsAa", "seAgVsAa")
write.csv(histEst, file = "historicalMultinomialModelEstimates3m.csv")

# save model fit
saveRDS(histMod, "histMultinomialGLMfit.rds")

# time-varying models ----------------------------------------------------------

# species frequency counts
sp1m <- cbind(dat1m$Aa, dat1m$Ac, dat1m$Ag)

# model selection covariates  --------------------------------------------------

seldat <- dat1m[ , c("lon.centroid", "lat.centroid",
                    "d2c", "d2r", "pop", "elev", "origin.year",
                    "PRECTOTCORR1m", "T2M1m", "RH2M1m")]

# use absolute value of latitude for insolation proxy
seldat$lat.centroid <- abs(seldat$lat.centroid)

## model including only 1 month met covariates ---------------------------------

# drop one month met covariates
m1Udat <- seldat

### untransformed covariates ---------------------------------------------------

# scale covariate data to speed up convergence
m1UdatScaled <- as.data.frame(apply(m1Udat, 2, function(x) (x - min(x))/(max(x) - min(x))))
raScaling.ls <- list(
  mins = apply(m1Udat, 2, min),
  maxs = apply(m1Udat, 2, max)
)
saveRDS(raScaling.ls, file = "raScaling.rds")

# fit initial full model
m1U.init <- multinom(sp1m ~ (lon.centroid + lat.centroid + d2r + pop + elev + d2c +
                               PRECTOTCORR1m + T2M1m + RH2M1m)^2 +
                       + I(lon.centroid^2) + I(lat.centroid^2) + I(d2r^2) + I(pop^2) + I(elev^2) + I(d2c^2) +
                       + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                       + origin.year, data = m1UdatScaled, maxit = 2000)

# backwards stepwise selection based on AIC
m1U.stp <- MASS::stepAIC(m1U.init,
                         scope = list(
                           upper = ~ (lon.centroid + lat.centroid + d2r + pop + elev + d2c +
                                        PRECTOTCORR1m + T2M1m + RH2M1m)^2 +
                             + I(lon.centroid^2) + I(lat.centroid^2) + I(d2r^2) + I(pop^2) + I(elev^2) + I(d2c^2) +
                             + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                             + origin.year,
                           lower = ~1),
                         direction = "backward", k = log(nrow(m1UdatScaled))
)

## transformed covariates ------------------------------------------------------

# log1p transform all nonnegative covariates (excluding longitude)
m1Tdat <- m1Udat
m1Tdat[ , -1] <- log1p(m1Udat[ , -1])

# scale covariate data to speed up convergence
m1TdatScaled <- as.data.frame(apply(m1Tdat, 2, function(x) (x - min(x))/(max(x) - min(x))))

# fit initial full model
m1T.init <- multinom(sp1m ~ (lon.centroid + lat.centroid + d2r + pop + elev + d2c +
                               PRECTOTCORR1m + T2M1m + RH2M1m)^2 +
                       + I(lon.centroid^2) + I(lat.centroid^2) + I(d2r^2) + I(pop^2) + I(elev^2) + I(d2c^2) +
                       + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                       + origin.year, data = m1TdatScaled, maxit = 2000)

# backwards stepwise selection based on AIC
m1T.stp <- MASS::stepAIC(m1T.init,
                         scope = list(
                           upper = ~ (lon.centroid + lat.centroid + d2r + pop + elev + d2c +
                                        PRECTOTCORR1m + T2M1m + RH2M1m)^2 +
                             + I(lon.centroid^2) + I(lat.centroid^2) + I(d2r^2) + I(pop^2) + I(elev^2) + I(d2c^2) +
                             + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                             + origin.year,
                           lower = ~1),
                         direction = "backward", k = log(nrow(m1TdatScaled))
)

## comparison ------------------------------------------------------------------

# BIC for initial full models
# best fitting model is with untransformed covariates
BIC(m1U.init, m1T.init)
BIC(m1U.stp, m1T.stp)

# selected final model ---------------------------------------------------------
# model with  untransformed covariates
ra.model <- multinom(sp1m ~ lon.centroid + lat.centroid + d2r + pop + elev + d2c +
                       PRECTOTCORR1m + T2M1m + RH2M1m + I(lon.centroid^2) + I(d2r^2) +
                       I(d2c^2) + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                       origin.year + lon.centroid:lat.centroid + lon.centroid:d2r +
                       lon.centroid:pop + lon.centroid:elev + lon.centroid:d2c +
                       lon.centroid:PRECTOTCORR1m + lon.centroid:T2M1m + lon.centroid:RH2M1m +
                       lat.centroid:d2r + lat.centroid:pop + lat.centroid:elev +
                       lat.centroid:d2c + lat.centroid:PRECTOTCORR1m + lat.centroid:T2M1m +
                       d2r:pop + d2r:elev + d2r:d2c + d2r:PRECTOTCORR1m + d2r:T2M1m +
                       d2r:RH2M1m + pop:elev + pop:d2c + pop:T2M1m + pop:RH2M1m +
                       elev:d2c + elev:T2M1m + elev:RH2M1m + d2c:PRECTOTCORR1m +
                       d2c:T2M1m + d2c:RH2M1m + PRECTOTCORR1m:T2M1m + PRECTOTCORR1m:RH2M1m +
                       T2M1m:RH2M1m,
                     data = m1UdatScaled, maxit = 1000)
all.equal(formula(m1U.stp), formula(ra.model))

# coefficients
raEst <- cbind(t(summary(ra.model)$coefficients), t(summary(ra.model)$standard.errors))
colnames(raEst) <- c("coeffAcVsAa", "coeffAgVsAa", "seAcVsAa", "seAgVsAa")
write.csv(raEst, file = "MultinomialModelEstimates1m.csv")

relABtab <- xtable::xtable(
  raEst,
  align = c("p{4.5cm}", "p{2.5cm}", "p{2cm}", "p{2.5cm}", "p{2cm}"),
  label = "label here",
  caption = "caption here",
  display = c("s", "f", "f", "f", "f")
)
cat(print(relABtab, math.style.exponents = TRUE),
    file = "relABtab.txt")


# save model object
saveRDS(ra.model, "multinomialGLMfit1m.rds")

# model fit
ra.pars <- list(
  coef = coef(ra.model),
  vcov = vcov(ra.model)
)
saveRDS(ra.pars, file = "raPars.rds")

## evaluation of model fit -----------------------------------------------------

### plot residuals -------------------------------------------------------------

# data.frame uses untransformed geographic coordinates
m.res <- as.data.frame(ra.model$residuals)
names(m.res) <- c("Aa", "Ac", "Ag")
m.res$long <- dat1m$lon.centroid
m.res$lat <- dat1m$lat.centroid

library(sf)
library(rnaturalearth)
library(ggplot2)

# sf object for residuals
res.sf <- st_as_sf(m.res,
                   coords = c("long", "lat"),
                   crs = "epsg:4326")

Africa <- ne_countries(continent = 'Africa', scale = 'medium', returnclass = 'sf')
Africa <- Africa[, c('name', 'name_long', 'iso_a3')]

pAa <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = res.sf, aes(colour = (Aa))) +
  scale_colour_distiller(type = "div") +
  xlim(-29, 50) + ylim(-30, 20)

pAc <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = res.sf, aes(colour = (Ac))) +
  scale_colour_distiller(type = "div") +
  xlim(-29, 50) + ylim(-30, 20)

pAg <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = res.sf, aes(colour = (Ag))) +
  scale_colour_distiller(type = "div") +
  xlim(-29, 50) + ylim(-30, 20)

png("residualsRA.png",
    width = 480*4, height = 480*2)
gridExtra::grid.arrange(pAa, pAc, pAg, ncol = 3)
dev.off()

### predictions versus observed ------------------------------------------------

mPred <- readRDS("multinomialGLMfit1m.rds")

# untransformed covariates for prediction
Preds.df <- m1Udat
# using max origin year
Preds.df$origin.year <- max(Preds.df$origin.year)
m.P <- predict(mPred, Preds.df, type = "probs")
colnames(m.P) <- c("Aa", "Ac", "Ag")

m.p.df <- as.data.frame(m.P)
# untransformed coordinates
m.p.df$long <- dat1m$lon.centroid
m.p.df$lat <- dat1m$lat.centroid

#### predictions ---------------------------------------------------------------

# quick plotting of predictions
# library(ggplot2)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(sf)

# sf object for probability predictions
m.p.sf <- st_as_sf(m.p.df,
                   coords = c("long", "lat"),
                   crs = "epsg:4326")

# bubbles by species predicted probs
predAa <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.p.sf, aes(size = (Aa)), pch = 20, colour = "red")

predAc <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.p.sf, aes(size = (Ac)), pch = 20, colour = "red")

predAg <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.p.sf, aes(size = (Ag)), pch = 20, colour = "red")

png("predictionsRA.png",
    width = 480*4, height = 480*2)
gridExtra::grid.arrange(predAa, predAc, predAg, ncol = 3)
dev.off()

#### observations -------------------------------------------------------------

odat <- cbind(dat1m[ , c("Aa", "Ac", "Ag")]/rowSums(dat1m[ , c("Aa", "Ac", "Ag")]),
              dat1m$lon.centroid, dat1m$lat.centroid)
names(odat)[4:5] <- c("long",  "lat")
m.o.sf <- st_as_sf(odat,
                   coords = c("long", "lat"),
                   crs = "epsg:4326")
# bubbles by species predicted probs
obsAa <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.o.sf, aes(size = (Aa)), pch = 20, colour = "red")

obsAc <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.o.sf, aes(size = (Ac)), pch = 20, colour = "red")

obsAg <- ggplot(Africa) +
  geom_sf() +
  geom_sf(data = m.o.sf, aes(size = (Ag)), pch = 20, colour = "red")

png("observationsRA.png",
    width = 480*4, height = 480*2)
gridExtra::grid.arrange(obsAa, obsAc, obsAg, ncol = 3)
dev.off()


