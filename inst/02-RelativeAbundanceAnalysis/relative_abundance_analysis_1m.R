# relative abundance analysis
# multinomial model
# one month only analysis

dat1m <- readRDS("../01-RelativeAbundanceData/ra1m_obs_covs_stra.rds")

library(nnet)

# time-varying models ----------------------------------------------------------

# species frequency counts
sp1m <- cbind(dat1m$Aa, dat1m$Ac, dat1m$Ag)

# model selection covariates  --------------------------------------------------

m1Udat <- dat1m[ , c("lon.centroid", "lat.centroid",
                    "d2c", "ppw", "elev", "origin.year", "hpd",
                    "PRECTOTCORR1m", "T2M1m", "RH2M1m", "IRS", "ITN")]

# use absolute value of latitude for insolation proxy
m1Udat$lat.centroid <- abs(m1Udat$lat.centroid)
# human population urban threshold
m1Udat$hpd <- ifelse(m1Udat$hpd < 2000, m1Udat$hpd, 2000)

# scale covariate data to speed up convergence ---------------------------------
m1UdatScaled <- as.data.frame(apply(m1Udat, 2, function(x) (x - min(x))/(max(x) - min(x))))
raScaling.ls <- list(
  mins = apply(m1Udat, 2, min),
  maxs = apply(m1Udat, 2, max)
)
saveRDS(raScaling.ls, file = "raScaling.rds")

# IRS intervention more likely regional and outside of An. coluzzii range
colSums(dat1m[ , c("Aa", "Ac", "Ag")])/sum((dat1m[ , c("Aa", "Ac", "Ag")]))
colSums(dat1m[dat1m$IRS > 0, c("Aa", "Ac", "Ag")])/sum(dat1m[dat1m$IRS > 0, c("Aa", "Ac", "Ag")])
# whereas ITN coverage occurs across all species ranges
colSums(dat1m[dat1m$ITN > 0, c("Aa", "Ac", "Ag")])/sum(dat1m[dat1m$ITN > 0, c("Aa", "Ac", "Ag")])

# fit initial full model with vector intervention spatio-temporal covariates
# more extensive ITN coverage: include main effect and year interaction
# less extensive IRS coverage: main IRS effect only because interaction may be
# conflated with taxonomic ID annual trend for Ac
m1U.init <- multinom(sp1m ~ (lon.centroid + lat.centroid + ppw + elev + d2c +
                               PRECTOTCORR1m + T2M1m + RH2M1m + hpd)^2 +
                       + I(lon.centroid^2) + I(lat.centroid^2) + I(ppw^2) + I(elev^2) + I(d2c^2) + I(hpd^2) + 
                       + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                       + (origin.year + ITN)^2,
                     data = m1UdatScaled, maxit = 2000)
BIC(m1U.init)

# backwards stepwise selection based on BIC
m1U.stp <- MASS::stepAIC(m1U.init,
                         scope = list(
                           upper = ~ (lon.centroid + lat.centroid + ppw + elev + d2c +
                                        PRECTOTCORR1m + T2M1m + RH2M1m + hpd)^2 +
                             + I(lon.centroid^2) + I(lat.centroid^2) + I(ppw^2) + I(elev^2) + I(d2c^2) + I(hpd^2) + 
                             + I(PRECTOTCORR1m^2) + I(T2M1m^2) + I(RH2M1m^2) +
                             + (origin.year + ITN)^2,
                           lower = ~ 1),
                         direction = "backward", k = log(sum(sp1m))
)

## comparison ------------------------------------------------------------------

# BIC for initial full models
BIC(m1U.init, m1U.stp)

# selected final model ---------------------------------------------------------
# model with  untransformed covariates
ra.model <- multinom(sp1m ~ lon.centroid + lat.centroid + ppw + elev + d2c + PRECTOTCORR1m + 
                       T2M1m + RH2M1m + hpd + I(lon.centroid^2) + I(lat.centroid^2) + 
                       I(ppw^2) + I(elev^2) + I(hpd^2) + I(PRECTOTCORR1m^2) + I(T2M1m^2) + 
                       I(RH2M1m^2) + origin.year + ITN + lon.centroid:lat.centroid + 
                       lon.centroid:elev + lon.centroid:d2c + lon.centroid:PRECTOTCORR1m + 
                       lon.centroid:T2M1m + lon.centroid:RH2M1m + lon.centroid:hpd + 
                       lat.centroid:ppw + lat.centroid:elev + lat.centroid:d2c + 
                       lat.centroid:T2M1m + lat.centroid:hpd + ppw:elev + ppw:d2c + 
                       ppw:PRECTOTCORR1m + ppw:T2M1m + ppw:RH2M1m + ppw:hpd + elev:d2c + 
                       elev:PRECTOTCORR1m + elev:T2M1m + elev:hpd + d2c:PRECTOTCORR1m + 
                       d2c:RH2M1m + d2c:hpd + PRECTOTCORR1m:T2M1m + PRECTOTCORR1m:RH2M1m + 
                       T2M1m:RH2M1m + T2M1m:hpd + RH2M1m:hpd + origin.year:ITN,
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
Preds.df <- m1UdatScaled
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


