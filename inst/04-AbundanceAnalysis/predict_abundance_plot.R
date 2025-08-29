# per capita abundance prediction plot

m <- readRDS("abundanceGLMfit.rds")
s.ls <- readRDS("abundanceScaling.rds")

# prediction design matrix
pop.pred <- seq(1, 70000, length = 101)
hpd.pred <- pop.pred/25 # (assuming 5 km x 5 km cell)
precip.pred <- seq(0, 30, length = 101)
X <- matrix(0, nrow = length(pop.pred), ncol = length(coef(m)),
                dimnames = list(
                  NULL, names(coef(m))
                ))
# populate base prediction design matrix at midpoints of covariates
X[ , "(Intercept)"] <- 1
# main effects
X[ , c("precip", "rh", "hpd", "itn")] <- 0.5
# squared terms and interactions
X[ , c("I(precip^2)", "I(hpd^2)", "precip:hpd", "rh:hpd")] <- 0.5^2

# prediction design matrix varying human population
hpd.cov <- ifelse(hpd.pred > 2000, 2000, hpd.pred) # threshold at 2K
Xpop <- X
Xpop[ , "hpd"] <- (hpd.cov - s.ls$mins["hpd"])/(s.ls$maxs["hpd"] - s.ls$mins["hpd"])
Xpop[ , "I(hpd^2)"] <- Xpop[ , "hpd"]^2
Xpop[ , "precip:hpd"] <- Xpop[ , "precip"]*Xpop[ , "hpd"]
Xpop[ , "rh:hpd"] <- Xpop[ , "rh"]*Xpop[ , "hpd"]
Xpop.df <- as.data.frame(Xpop)
Xpop.df$units <- 1

# prediction design matrix varying precipitation
Xprecip <- X
Xprecip[ , "precip"] <- (precip.pred - s.ls$mins["precip"])/(s.ls$maxs["precip"] - s.ls$mins["precip"])
Xprecip[ , "I(precip^2)"] <- Xprecip[ , "precip"]^2
Xprecip[ , "precip:hpd"] <- Xprecip[ , "precip"]*Xprecip[ , "hpd"]
Xprecip.df <- as.data.frame(Xprecip)
Xprecip.df$units <- 1

# predictions
hpd.predict <- predict(m, newdata = Xpop.df, type = "response", se.fit = TRUE)
Xtot.df <- Xpop.df
Xtot.df$units <- ifelse(pop.pred > 50000, 50000, pop.pred) # theshold at 50K
tot.predict <- predict(m, newdata = Xtot.df, type = "response", se.fit = TRUE)
precip.predict <- predict(m, newdata = Xprecip.df, type = "response", se.fit = TRUE)

# plots 
options(scipen=10000)
png("abundance.png", width = 480*2, height = 480*3, pointsize = 12*3)

par(mfrow = c(3, 1), mar = c(5, 6, 1, 1), mgp = c(4.5, 1, 0))

# population (per capita)
plot(hpd.pred, hpd.predict$fit, type = 'l', ylim = c(0, 45), las = 1,
     xlab = "",
     ylab = expression(paste("Adult female ", italic(An.~gambiae~s.l.))))
lines(hpd.pred, hpd.predict$fit + hpd.predict$se.fit, lty = 2)
lines(hpd.pred, hpd.predict$fit - hpd.predict$se.fit, lty = 2)
mtext("Human population per square kilometre", side = 1, line = 3, cex = 0.8)
mtext("(per capita)", side = 2, line = 3.25, cex = 2/3)
abline(v = 2000, lty = 3)

# population (total)
plot(pop.pred, tot.predict$fit, type = 'l', ylim = c(0, 6e5), las = 1,
     xlab = "",
     ylab = expression(paste("Adult female ", italic(An.~gambiae~s.l.))))
lines(pop.pred, tot.predict$fit + tot.predict$se.fit, lty = 2)
lines(pop.pred, tot.predict$fit - tot.predict$se.fit, lty = 2)
mtext("Human population per 5 km x 5 km cell", side = 1, line = 3, cex = 0.8)
abline(v = 50000, lty = 3)

# precipitation (per capita)
plot(precip.pred, precip.predict$fit, type = 'l', ylim = c(0, 20), las = 1,
     xlab = "",
     ylab = expression(paste("Adult female ", italic(An.~gambiae~s.l.))))
lines(precip.pred, precip.predict$fit + precip.predict$se.fit, lty = 2)
lines(precip.pred, precip.predict$fit - precip.predict$se.fit, lty = 2)
mtext("Daily precipitation (mm)", side = 1, line = 3, cex = 0.8)
mtext("(per capita)", side = 2, line = 3.25, cex = 2/3)

dev.off()

# midpoints
0.5*(s.ls$maxs["hpd"] - s.ls$mins["hpd"]) + s.ls$mins["hpd"]
0.5*(s.ls$maxs["precip"] - s.ls$mins["precip"]) + s.ls$mins["precip"]
0.5*(s.ls$maxs["rh"] - s.ls$mins["rh"]) + s.ls$mins["rh"]
0.5*(s.ls$maxs["itn"] - s.ls$mins["itn"]) + s.ls$mins["itn"]

