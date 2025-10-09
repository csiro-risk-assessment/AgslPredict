# xi estimation

# abundance data ----------------------------------------------

abund <- readRDS("../03-AbundanceData/abund_obs_covs.rds")

# original data for model
dat <- abund

# observations and effort
obsMat <- dat[ , c("Agsl", "units")]

# collection method submodel ----------------------------------------------

# collection codes
# combine animal catches
# light traps catches grouped except for where denoted as outdoor
codes.df <- data.frame(
  protocols = c(
    "animal baited net trap catch",       "animal biting catch - outdoors",     "artificial pit shelter",
    "CDC light trap catch",               "CDC miniature light trap model 512",
    "exit trap catch",                    "indoor light trap catch",            "man biting catch",
    "man biting catch - indoors",         "man biting catch - outdoors",        "outdoor light trap catch",
    "pyrethrum spray catch"
  ),
  codes = c("ABC", "ABC",  "APS",  "LTC",  "LTC", "ETC",
            "ILTC",  "MBC",   "MBCI",  "MBCO",  "OLTC",  "PSC")
)

# coding matrix for collection methods
# treatment contrasts
# known zero intercept at indoor HLC
collectionCodes <- sort(unique(codes.df$codes))
# design matrix for collection
collectMat <- matrix(0, nrow = length(collectionCodes), ncol = 6,
                     dimnames = list(
                       code = collectionCodes,
                       factor = c("light", "exit", "PSC", "animal", "pit",
                                  "outdoor")
                     )
)
# light
collectMat[grepl("LTC", collectionCodes), "light"] <- 1
# exit
collectMat[grepl("ETC", collectionCodes), "exit"] <- 1
# PSC
collectMat[grepl("PSC", collectionCodes), "PSC"] <- 1
# animal
collectMat[grepl("ABC", collectionCodes), "animal"] <- 1
# pit
collectMat[grepl("APS", collectionCodes), "pit"] <- 1
# outdoor unknown - will overwrite MBC and LTC
collectMat[grepl("MBC", collectionCodes), "outdoor"] <-
  collectMat[grepl("LTC", collectionCodes), "outdoor"] <- 0.5
# outdoor known
collectMat[grepl("ABC", collectionCodes), "outdoor"] <-
  collectMat[grepl("APS", collectionCodes), "outdoor"] <-
  collectMat[grepl("MBCO", collectionCodes), "outdoor"] <-
  collectMat[grepl("OLTC", collectionCodes), "outdoor"] <- 1
# indoor known
collectMat[grepl("MBCI", collectionCodes), "outdoor"] <-
  collectMat[grepl("ILTC", collectionCodes), "outdoor"] <- 0

# design matrices for collection methods, AcAg vs Aa
collectMatAgsl <- matrix(NA, nrow = nrow(dat), ncol = ncol(collectMat))
colnames(collectMatAgsl) <- colnames(collectMat)
rownames(collectMatAgsl) <- dat$Collection.protocols
for (i in 1:length(collectionCodes)) {
  u.collections <- codes.df[codes.df$codes == collectionCodes[i], "protocols"]
  for (jj in 1:length(u.collections)) {
    cm <- matrix(collectMat[collectionCodes[i], ],
                 nrow = sum(rownames(collectMatAgsl) == u.collections[jj]),
                 ncol = ncol(collectMatAgsl), byrow = TRUE)
    collectMatAgsl[rownames(collectMatAgsl) == u.collections[jj], ] <- cm
  }
}

# q submodel ----------------------------------------------

# design matrix
# population density covariate thresholded at 2000/km^2
#   - 2K/km^2 --> at least 50000 total population in 25 km^2 cell
#   - total population over 50K in contiguous 1 km^2 cells of greater than 
#     1500/km^2 density is high density urban area (World Bank 2020)
# elev
# abs lat
# daily precip, relative humidity
# ITN
qMat <- cbind(
  ifelse(dat$hpd < 2000, dat$hpd, 2000), 
  dat$elev,
  abs(dat$lat.centroid),
  rowMeans(dat[ , grepl("p_lag", names(dat))]),
  rowMeans(dat[ , grepl("r_lag", names(dat))]),
  dat$ITN
)
colnames(qMat) <- c(
  "hpd", "elev", "lat", 
  "precip", "rh",
  "itn"
)

# unit scalings ---------------------------------------------

sMat <- apply(qMat, 2, function(x) (x - min(x))/(max(x) - min(x)))
scaling.ls <- list(
  mins = apply(qMat, 2, min),
  maxs = apply(qMat, 2, max)
)
saveRDS(scaling.ls, file = "abundanceScaling.rds")

# analysis ---------------------------------------------------------------------

X <- cbind(
  as.matrix(collectMatAgsl),
  sMat
)
colnames(X) <- c(
  colnames(collectMatAgsl),
  colnames(sMat)
)

df <- as.data.frame(cbind(
  obsMat[ , c("Agsl", "units")],
  X
))

# Negative binomial with optimised theta
library(glm2)

fit.init <- MASS::glm.nb(Agsl ~ light + exit + PSC + animal + pit + outdoor +
                            precip + I(precip^2) +
                            rh + 
                            precip:rh +
                            precip:hpd +
                            hpd:rh +
                            hpd + I(hpd^2) +
                            itn +
                            offset(log(units)),
                          method = "glm.fit2",
                          maxit = 1000,
                          data = df)
fit.stp <- MASS::stepAIC(fit.init,
                         scope = list(
                           upper = ~ light + exit + PSC + animal + pit + outdoor +
                             precip + I(precip^2) +
                             rh + 
                             precip:rh +
                             precip:hpd +
                             hpd:rh +
                             hpd + I(hpd^2) +
                             itn +
                             offset(log(units)),
                           lower = ~ light + exit + PSC + animal + pit + outdoor +
                             precip + I(precip^2) +
                             rh + 
                             hpd + I(hpd^2) +
                             itn +
                             offset(log(units))),
                         direction = "backward", k = log(nrow(df))
)
BIC(fit.init, fit.stp)

fit.model <- MASS::glm.nb(Agsl ~ light + exit + PSC + animal + pit + outdoor + precip + 
                                    I(precip^2) + rh + hpd + I(hpd^2) + itn + precip:hpd + rh:hpd + 
                                    offset(log(units)),                           
                                    method = "glm.fit2",
                                    maxit = 1000,
                                    data = df)
all.equal(formula(fit.model), formula(fit.stp))
fit.model$converged
summary(fit.model)
BIC(fit.model)

theta.est <- fit.model$theta

## overdispersed model summary -------------------------------------------------

# coefficients
summary(fit.model)

# save model object
saveRDS(fit.model, "abundanceGLMfit.rds")

library(texreg)

texreg(fit.model, file = "negbinom_table.txt",
       caption = paste0("Negative binomial GLM with estimated overdipersion parameter of ",
                        round(theta.est, 3), "."))

