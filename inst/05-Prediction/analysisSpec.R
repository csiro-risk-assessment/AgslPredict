# read in thresholded (non urban) predictions ----------------------------------
# using results from prediction_plots.R
results.nu <- readRDS("resultsNU.rds")

# grand mean
xbar <- colMeans(results.nu[ , c("Aa_mean", "Ac_mean", "Ag_mean")], na.rm = TRUE)
names(xbar) <- c("Aa", "Ac", "Ag")

# dropping zero abundance sites
x <- results.nu[ , c("Aa_mean", "Ac_mean", "Ag_mean")]
x <- x[!is.na(x[ , 1]), ]
x <- x[x[ , 1] > 0, ] # all positive: conditional on non-zero abundances
x.bar <- colMeans(x)
p.bar <- x.bar/sum(x.bar)

# Pianka
r <- x/rowSums(x)
P <- diag(3)
P[1, 2] <- sum(r[ , 1]*r[ , 2])/sqrt(sum(r[ , 1]^2)*sum(r[ , 2]^2))
P[1, 3] <- sum(r[ , 1]*r[ , 3])/sqrt(sum(r[ , 1]^2)*sum(r[ , 3]^2))
P[2, 3] <- sum(r[ , 2]*r[ , 3])/sqrt(sum(r[ , 2]^2)*sum(r[ , 3]^2))
P[lower.tri(P)] <- P[upper.tri(P)]

# density independent parameters -----------------------------------------------

# Table 3.3, Hosack et al. 2023
sigma <- c(0.85, 0.84, 0.84)
names(sigma) <- c("Aa", "Ac", "Ag")

# Table C.1, Hosack et al. 2023
lambda <- 9
z <- 1 - sigma
b <- (lambda - z)/z

# Jacobian ---------------------------------------------------------------------
Jfun <- function(xbar, A, delta = 10) {

  x <- xbar

  # build Jacobian matrix
  n.age <- delta + 1 # include age zero
  J <- matrix(0, nrow = 3*(delta + 1), ncol = 3*(delta + 1))

  # lags
  J[2:n.age, 1:(n.age - 1)] <-
    J[(n.age + 2):(2*n.age), (n.age + 1):(2*n.age - 1)] <-
    J[(2*n.age + 2):(3*n.age), (2*n.age + 1):(3*n.age - 1)] <-
    diag(delta)

  # density independent adult mortality
  J[n.age, n.age] <- sigma[1]
  J[2*n.age, 2*n.age] <- sigma[2]
  J[3*n.age, 3*n.age] <- sigma[3]

  # density dependence
  # rowwise
  # species 1
  J[1, n.age] <- z[1]*(1 - z[1]*x[1]/(lambda*q[1]))
  J[1, 2*n.age] <- -z[1]^2*A[1, 2]*x[1]/(lambda*q[1])
  J[1, 3*n.age] <- -z[1]^2*A[1, 3]*x[1]/(lambda*q[1])
  # species 2
  J[n.age + 1, n.age] <- -z[2]^2*A[2, 1]*x[2]/(lambda*q[2])
  J[n.age + 1, 2*n.age] <- z[2]*(1 - z[2]*A[2, 2]*x[2]/(lambda*q[2]))
  J[n.age + 1, 3*n.age] <- -z[2]^2*A[2, 3]*x[2]/(lambda*q[2])
  # species 3
  J[2*n.age + 1, n.age] <- -z[3]^2*A[3, 1]*x[3]/(lambda*q[3])
  J[2*n.age + 1, 2*n.age] <- -z[3]^2*A[3, 2]*x[3]/(lambda*q[3])
  J[2*n.age + 1, 3*n.age] <- z[3]*(1 - z[3]*A[3, 3]*x[3]/(lambda*q[3]))

  J

}

## simulation ------------------------------------------------------------------

# Pianka
A <- P
q <- diag(z)%*%solve(diag(lambda - z))%*%A%*%x.bar
solve(A)%*%diag(lambda - z)%*%solve(diag(z))%*%q - x.bar

J <- Jfun(x.bar, A)
summary(abs(eigen(J)$values))

set.seed(999)
# initial condition: random deviation from overall average
xi0 <- rlnorm(3, log(x.bar), sdlog = 1)

Q <- diag(c(q))
N <- 150
xi <- matrix(NA, nrow = 3, ncol = N)
xi[ , 1:10] <- matrix(xi0, nrow = 10, ncol = 3)
for (tt in 11:N) {
  xi[ , tt] <- diag(sigma)%*%xi[ , tt - 1] + solve(diag(c(1 + solve(Q)%*%A%*%xi[ , tt - 10])))%*%xi[ , tt - 10]*lambda
}

plot(c(1, ncol(xi)), range(xi), ylog = TRUE, type = 'n',
     xlab = "Time", ylab = "Abundance")
lines(xi[1, ])
lines(xi[2, ], lty = 2)
lines(xi[3, ], lty = 3)
points(rep(ncol(xi), 3), x.bar) # overall average
legend('bottomright', bty = 'n',
       lty = 1:3, legend = c("Aa", "Ac", "Ag"))

## plot basin  -----------------------------------------------------------------

Jpfun <- function(p2, p3, A, delta = 10) {

  # proportions
  if (p2 < 0 || p3 < 0) stop("req positive proportions")
  if (p2 > 1 || p3 > 1) stop("req proportions lt one")

  p1 <- 1 - p2 - p3

  x <- c(1, p2/p1, p3/p1)

  # build Jacobian matrix
  n.age <- delta + 1 # include age zero
  J <- matrix(0, nrow = 3*(delta + 1), ncol = 3*(delta + 1))

  # lags
  J[2:n.age, 1:(n.age - 1)] <-
    J[(n.age + 2):(2*n.age), (n.age + 1):(2*n.age - 1)] <-
    J[(2*n.age + 2):(3*n.age), (2*n.age + 1):(3*n.age - 1)] <-
    diag(delta)

  # density independent adult mortality
  J[n.age, n.age] <- sigma[1]
  J[2*n.age, 2*n.age] <- sigma[2]
  J[3*n.age, 3*n.age] <- sigma[3]

  # density dependence
  # rowwise
  # species 1
  J[1, n.age] <- z[1]*(1 - ((lambda - z[1])/lambda)*x[1]/(A[1, ]%*%x))
  J[1, 2*n.age] <- -z[1]*((lambda - z[1])/lambda)*A[1, 2]*x[1]/(A[1, ]%*%x)
  J[1, 3*n.age] <- -z[1]*((lambda - z[1])/lambda)*A[1, 3]*x[1]/(A[1, ]%*%x)
  # species 2
  J[n.age + 1, n.age] <- -z[2]*((lambda - z[2])/lambda)*A[2, 1]*x[2]/(A[2, ]%*%x)
  J[n.age + 1, 2*n.age] <- z[2]*(1 - ((lambda - z[2])/lambda)*x[2]/(A[2, ]%*%x))
  J[n.age + 1, 3*n.age] <- -z[2]*((lambda - z[2])/lambda)*A[2, 3]*x[2]/(A[2, ]%*%x)
  # species 3
  J[2*n.age + 1, n.age] <- -z[3]*((lambda - z[3])/lambda)*A[3, 1]*x[3]/(A[3, ]%*%x)
  J[2*n.age + 1, 2*n.age] <- -z[3]*((lambda - z[3])/lambda)*A[3, 2]*x[3]/(A[3, ]%*%x)
  J[2*n.age + 1, 3*n.age] <- z[3]*(1 - ((lambda - z[3])/lambda)*x[3]/(A[3, ]%*%x))

  J

}

## plot ------------------------------------------------------------------------

# response functions use global vars
FuncP <- function(a, b, c) {
  Jint <- Jpfun(b, c, P)
  abs.eigs <- abs(eigen(Jint)$values)
  max(abs.eigs)
}

pts <- expand.grid(a = seq(0.01, 0.99, len = 101),
                   b = seq(0.01, 0.99, len = 101))
pts <- pts[rowSums(pts) < 1, ]
pts$c <- 1 - rowSums(pts)

# plots use Ternary package
library(Ternary)

tri <- TriangleCentres(resolution = 50L)
tern.xy <- XYToTernary(tri["x", ], tri["y", ])

vP <- numeric(ncol(tern.xy))
for (i in 1:ncol(tern.xy)) {
  vP[i] <- FuncP(tern.xy["a", i], tern.xy["b", i], tern.xy["c", i])
}

png("srtPianka.png", width = 480*1, height = 480,
    pointsize = 24)
par(mar = rep(0.2, 4))

TernaryPlot(alab = expression(paste("Percentage of  ", italic(An.~arabiensis))),
            blab = expression(paste("Percentage of  ", italic(An.~coluzzii))),
            clab = expression(paste("Percentage of  ", italic(An.~gambiae), "  s.s.")))
mapP <- rbind(x = tri["x", ], y = tri["y", ], z = vP,
                down = tri["triDown", ])
ColourTernary(mapP, spectrum = viridisLite::viridis(256L, alpha = 0.6))

PlotTools::SpectrumLegend(
  "topleft",
  legend = rev(round(seq(min(vP), max(vP), length.out = 4), 3)),
  palette = viridisLite::viridis(256L, alpha = 0.6),
  bty = "n",    # No framing box
  inset = 0.02,
  xpd = NA      # Do not clip at edge of figure
)

# estimate avg equilibrium
TernaryPoints(p.bar, pch = 4, cex = 1.7, col = "white")
# estimate min spec rad
TernaryPoints(tern.xy[1:3, which.min(vP)], pch = 1, cex = 1.7, col = "white")

dev.off()


# optimised analysis------------------------------------------------------------

obj <- function(p23, Ap) {
  Jp <- Jpfun(p23[1], p23[2], Ap)
  max(abs(eigen(Jp)$values))
}

objJ <- function(log.aij, p.bar) {

  aij <- exp(log.aij)
  Amat <- diag(3)
  Amat[lower.tri(Amat)] <- aij
  Amat[upper.tri(Amat)] <- Amat[lower.tri(Amat)]

  rhobar <- obj(p.bar[2:3], Amat)

  if (rhobar < 1) {

    res <- (rhobar - 1)*mean(exp(log.aij))

  } else {

    res <- Inf

  }

  res

}

# grid search
grid.pts <- expand.grid(a = seq(0, 1, len = 101),
                        b = seq(0, 1, len = 101),
                        c = seq(0, 1, len = 101))
grid.pts <- as.matrix(grid.pts)
grid.pts <- cbind(grid.pts, NA)
colnames(grid.pts)[ncol(grid.pts)] <- "res"
summary(grid.pts)

pt <- proc.time()
for (i in 1:nrow(grid.pts)) {
  grid.pts[i, "res"] <- objJ(log(c(grid.pts[i, 1:3])), p.bar)
  if (i%%1000 == 0) cat('completed ', 100*(i/nrow(grid.pts)), "percent after ", proc.time()[3] - pt[3], "secs \n")
}

ind.maxGrad <- which.min(grid.pts[ , "res"])
grid.pts[ind.maxGrad, ]
maxGrad.coords <- grid.pts[ind.maxGrad, 1:3]

Agrid <- diag(3)
Agrid[lower.tri(Agrid)] <- maxGrad.coords
Agrid[upper.tri(Agrid)] <- Agrid[lower.tri(Agrid)]
round(Agrid, 3)

print(xtable::xtable(Agrid), file = "Agrid.txt")

Aoptim <- Agrid

J.o.equal <- Jpfun(1/3, 1/3, Aoptim)
J.o.xbar <- Jpfun(p.bar[2], p.bar[3], Aoptim)

summary(abs(eigen(J.o.equal)$values))
summary(abs(eigen(J.o.xbar)$values))

## simulation ------------------------------------------------------------------

# optimised
A <- Aoptim
q <- diag(z)%*%solve(diag(lambda - z))%*%A%*%x.bar
solve(A)%*%diag(lambda - z)%*%solve(diag(z))%*%q - x.bar

J <- Jfun(x.bar, A)
all.equal(Jfun(x.bar, A), Jpfun(x.bar[2]/sum(x.bar), x.bar[3]/sum(x.bar), A))
summary(abs(eigen(J)$values))

set.seed(999)
# initial condition: random deviation from overall average
xi0 <- rlnorm(3, log(x.bar), sdlog = 1)

Q <- diag(c(q))
N <- 150*2
xi <- matrix(NA, nrow = 3, ncol = N)
xi[ , 1:10] <- matrix(xi0, nrow = 10, ncol = 3)
for (tt in 11:N) {
  xi[ , tt] <- diag(sigma)%*%xi[ , tt - 1] + solve(diag(c(1 + solve(Q)%*%A%*%xi[ , tt - 10])))%*%xi[ , tt - 10]*lambda
}

plot(c(1, ncol(xi)), range(xi), ylog = TRUE, type = 'n',
     xlab = "Time", ylab = "Abundance")
lines(xi[1, ])
lines(xi[2, ], lty = 2)
lines(xi[3, ], lty = 3)
points(rep(ncol(xi), 3), x.bar) # overall average
legend('bottomright', bty = 'n',
       lty = 1:3, legend = c("Aa", "Ac", "Ag"))


## plot ------------------------------------------------------------------------

# response functions use global vars
FuncOpt <- function(a, b, c) {
  Jint <- Jpfun(b, c, Aoptim)
  abs.eigs <- abs(eigen(Jint)$values)
  max(abs.eigs)
}

pts <- expand.grid(a = seq(0.01, 0.99, len = 101),
                   b = seq(0.01, 0.99, len = 101))
pts <- pts[rowSums(pts) < 1, ]
pts$c <- 1 - rowSums(pts)

# plots use Ternary package
library(Ternary)

tri <- TriangleCentres(resolution = 50L)
tern.xy <- XYToTernary(tri["x", ], tri["y", ])

vOpt <- numeric(ncol(tern.xy))
for (i in 1:ncol(tern.xy)) {
  vOpt[i] <- FuncOpt(tern.xy["a", i], tern.xy["b", i], tern.xy["c", i])
}

png("specOptim.png", width = 480*1, height = 480,
    pointsize = 24)
par(mar = rep(0.2, 4))

TernaryPlot(alab = expression(paste("Percentage of  ", italic(An.~arabiensis))),
            blab = expression(paste("Percentage of  ", italic(An.~coluzzii))),
            clab = expression(paste("Percentage of  ", italic(An.~gambiae), "  s.s.")))
mapOpt <- rbind(x = tri["x", ], y = tri["y", ], z = (vOpt),
             down = tri["triDown", ])
ColourTernary(mapOpt, spectrum = viridisLite::viridis(256L, alpha = 0.6))

PlotTools::SpectrumLegend(
  "topleft",
  legend = rev(round(seq(min(vOpt), max(vOpt), length.out = 4), 3)),
  palette = viridisLite::viridis(256L, alpha = 0.6),
  bty = "n",    # No framing box
  inset = 0.02,
  xpd = NA      # Do not clip at edge of figure
)

# estimate avg equilibrium
TernaryPoints(p.bar, pch = 4, cex = 1.7, col = "white")
# estimate min spec rad
TernaryPoints(tern.xy[1:3, which.min(vOpt)], pch = 1, cex = 1.7, col = "white")

dev.off()

# table ------------------------------------------------------------------------

Aij <- matrix(NA, nrow = 2, ncol = 3,
              dimnames = list(
                c("Pianka", "Optimised"),
                c("AaAc", "AaAg", "AcAg")
              ))
Aij[1, ] <- P[lower.tri(P)]
Aij[2, ] <- Aoptim[lower.tri(Aoptim)]
Aij <- round(Aij, 3)

library(xtable)
print(xtable(Aij), file = "Aij.txt")


