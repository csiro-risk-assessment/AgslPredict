# raPreds: data.frame of covariates for multinomial model
# aPreds: data.frame of covariates for abundance model
# m: multinomial GLM model fit for relative abundance
# ra.pars: multinomial parameter estimates and vcov
# a.mod: abundance GLM model fit
pxPrediction <- function(raPreds, aPreds, m, ra.pars, a.mod) {

  require(nnet)

  # data checks ----------------------------------------------------------------

  if (!all(class(m) == c("multinom", "nnet")))
    stop("expecting multinomial GLM model object")

  # data transformations for abundance predictions -----------------------------

  # multiplier for spatio-temporal total human population in cell
  human.pop <- aPreds$pop

  # total obs abundance prediction ---------------------------------------------

  xi <- predict(a.mod, aPreds, type = "response")
  # average prediction variance of linear predictor
  xi.v <- mean(predict(a.mod, aPreds, type = "link", se.fit = TRUE)$se.fit^2)

  # species composition prediction (observable rel abundance) ------------------
  r <- predict(m, raPreds, type = "probs")
  colnames(r) <- c("Aa", "Ac", "Ag")
  # average prediction variance of linear predictors
  X <- model.matrix(formula(m)[-2], raPreds)
  X1 <- cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(X)))
  X2 <- cbind(matrix(0, nrow = nrow(X), ncol = ncol(X)), X)
  r.v <- c(
      sum(diag(crossprod(X1, X1)%*%ra.pars$vcov)),
      sum(diag(crossprod(X2, X2)%*%ra.pars$vcov))
    )/nrow(X)
  names(r.v) <- c("logAcAa", "logAgAa")

  # proportional to actual rel abundance: --------------------------------------
  # r adjusted for human feeding rate
  xi.adj.hfr <-  xi/rowSums(r*matrix(c(0.37, 0.49, 0.47),
                            nrow = nrow(r), ncol = 3,
                            byrow = TRUE))

  # actual tot abundance adjust for human pop ----------------------------------
  x <- r*xi.adj.hfr*human.pop

  list(x = x, xi.m = mean(xi), xi.v = xi.v,
       r = colMeans(r), r.v = r.v)

}
