# xi estimation

# abundance data ----------------------------------------------

abund <- readRDS("../03-AbundanceData/abund_obs_covs.rds")

# original data for model
dat <- abund

# collection method management
# exclude larvae collections - targeting adult females
dat <- dat[dat$Collection.protocols != "collection of larvae from dippers", ]
# drop all records with collection methods having less than 10 occurrences
rareMethods <- names(table(dat$Collection.protocols))[table(dat$Collection.protocols) <= 10]
dat <- dat[!dat$Collection.protocols%in%rareMethods, ]

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
# log1p population
# log1p d2r
# log1p elev
# abs lat
# long
# daily precip, temp, RH
# collection method

qMat <- cbind(
  dat$d2c,
  dat$d2r,
  dat$pop,
  dat$elev,
  abs(dat$lat.centroid),
  dat$lon.centroid,
  rowMeans(dat[ , grepl("p_lag", names(dat))]),
  rowMeans(dat[ , grepl("t_lag", names(dat))]),
  rowMeans(dat[ , grepl("r_lag", names(dat))])
)
colnames(qMat) <- c(
  "d2c", "d2r", "pop", "elev", "lat", "long",
  "precip", "temp", "rh"
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
  "d2c", "d2r", "pop", "elev", "lat", "long",
  "precip", "temp", "rh"
)

df <- as.data.frame(cbind(
  obsMat[ , c("Agsl", "units")],
  X
))

# linear predictor
# precipitation: quadratic, interactions with all non obs process vars
# d2r: interaction with terrain covars: d2c and elev

# Poisson
fit.p <- glm(Agsl ~ light + exit + PSC + animal + pit + outdoor +
                       precip + I(precip^2) +
                       precip:lat + precip:elev +
                       precip:temp +
                       pop +
                       offset(log(units)),
                      family = poisson,
                      # method = "glm.fit2", # glm2 overwrite
                      maxit = 100,
                      data = df)
fit.p$converged

# Negative binomial with optimised theta
library(glm2)

# pop
fit.model <- MASS::glm.nb(Agsl ~ light + exit + PSC + animal + pit + outdoor +
                            precip + I(precip^2) +
                            precip:lat + precip:elev +
                            precip:temp +
                            pop +
                            offset(log(units)),
                          method = "glm.fit2",
                          data = df)
fit.model$converged
summary(fit.model)
curve(exp(coef(fit.model)[1] + coef(fit.model)["precip"]*x + coef(fit.model)["I(precip^2)"]*x^2)*x*1500)
curve(exp(coef(fit.model)[1] + coef(fit.model)["pop"]*x)*x*1500)

## evaluation of overdispersed model fit ---------------------------------------

BIC(fit.p, fit.model) # fitted theta = 0.232 neg bin substantially improves BIC
# est theta slightly below Bayesian estimate in Hosack et al. (2003)
# that uses different linear predictor (posterior mean: 0.26)
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

# historical -------------------------------------------------------------------

# q submodel ----------------------------------------------

# design matrix
# log1p population
# log1p d2r
# log1p elev
# abs lat
# long
# quadratic daily precip
# collection method

qMat <- cbind(
  log1p(dat$pop),
  log1p(dat$d2r),
  log1p(dat$elev),
  abs(dat$lat.centroid),
  dat$lon.centroid,
  dat[ , grepl("lag", names(dat))],
  dat[ , grepl("lag", names(dat))]^2
)
names(qMat) <- c(
  names(qMat)[!grepl("lag", names(qMat))],
  paste0("plag", 1:10),
  paste0("tlag", 1:10),
  paste0("rlag", 1:10),
  paste0("plagSq", 1:10),
  paste0("tlagSq", 1:10),
  paste0("rlagSq", 1:10)
)
# spatial covariates
qSpat <- qMat[ , !grepl("lag", names(qMat))]
# linear lagged weekly running mean precip
qPrecip <- qMat[ , grepl("plag", names(qMat)) & !grepl("Sq", names(qMat))]
# squared lagged weekly running mean precip
qPrecipSq <- qMat[ , grepl("plag", names(qMat)) & grepl("Sq", names(qMat))]

# export scalings ---------------------------------------------

scaling.ls <- list(
  spatial = list(
    centre = colMeans(qSpat),
    scaling = apply(qSpat, 2, sd)
  ),
  precip = list(
    centre = colMeans(qPrecip),
    scaling = apply(qPrecip, 2, sd)
  ),
  precipSq = list(
    centre = colMeans(qPrecipSq),
    scaling = apply(qPrecipSq, 2, sd)
  )
)

all.equal(scale(qSpat),
          sweep(as.matrix(qSpat), 2, scaling.ls$spatial$centre)%*%diag(1/scaling.ls$spatial$scaling),
          check.attributes = FALSE
)
all.equal(scale(qPrecip),
          sweep(as.matrix(qPrecip), 2, scaling.ls$precip$centre)%*%diag(1/scaling.ls$precip$scaling),
          check.attributes = FALSE
)
all.equal(scale(qPrecipSq),
          sweep(as.matrix(qPrecipSq), 2, scaling.ls$precipSq$centre)%*%diag(1/scaling.ls$precipSq$scaling),
          check.attributes = FALSE
)

sSpat <- sweep(as.matrix(qSpat), 2, scaling.ls$spatial$centre)%*%diag(1/scaling.ls$spatial$scaling)
sPrecip <- sweep(as.matrix(qPrecip), 2, scaling.ls$precip$centre)%*%diag(1/scaling.ls$precip$scaling)
sPrecipSq <- sweep(as.matrix(qPrecipSq), 2, scaling.ls$precipSq$centre)%*%diag(1/scaling.ls$precipSq$scaling)

# analysis ---------------------------------------------------------------------

X <- cbind(
  as.matrix(collectMatAgsl),
  sSpat,
  rowMeans(sPrecip),
  rowMeans(sPrecipSq)
)
colnames(X) <- c(
  colnames(collectMatAgsl),
  "pop", "d2r", "elev", "lat", "long",
  "precip", "precip.sq"
)

df <- as.data.frame(cbind(
  obsMat[ , c("Agsl", "units")],
  X
))

# Poisson model
fit.p <- glm(Agsl ~ -1 + light + exit + PSC + animal + pit + outdoor +
                          pop + d2r + elev + lat + long +
                          precip + precip.sq +
                          offset(log(units)),
                        family = "poisson",
                        method = "glm.fit2", # glm2 overwrite
                        maxit = 2000,
                        data = df)

# very close to Bayesian model in report (Hosack et al. 2023)
# that uses similar precip and elevation
fit.g22 <- MASS::glm.nb(Agsl ~ -1 + light + exit + PSC + animal + pit + outdoor +
                          pop + d2r + elev + lat + long +
                          precip + precip.sq +
                          offset(log(units)),
                        method = "glm.fit2", # glm2 overwrite
                        maxit = 2000,
                        data = df)
fit.g22$converged
summary(fit.g22)
# above Poisson model much worse fit compared to historical
# above neg bin model not quite as good fit compared to historical
BIC(fit.p, fit.g22)

saveRDS(fit.g22, "abundanceHistorical.rds")
