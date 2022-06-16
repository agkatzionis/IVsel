
#####################   IVSEL FUNCTIONS   #####################

## This contains the R code for the R fuctions we 
## need to run the "Instruments for selection" methods.

## Auxiliary expit and logit functions.
expit <- function (x) exp(x) / (1 + exp(x))
logit <- function (x) log(x / (1 - x))
expit2 <- function (x) {
  res <- exp(x) / (1 + exp(x))
  res[x > 0] <- 1 / (1 + exp(- x[x > 0]))
  res
}

## Inverse probability weighting - linear case.
ipw.linear <- function (X, Y, R, Z = NA, drop.x = FALSE) {
  
  ## Regress participation on covariates, compute scores.
  if (all(is.na(Z))) {
    ipw.par <- unname(glm(R ~ X, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X) %*% ipw.par))
  } else if (drop.x == TRUE) {
    ipw.par <- unname(glm(R ~ Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, Z) %*% ipw.par))
  } else {
    ipw.par <- unname(glm(R ~ X + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X, Z) %*% ipw.par))
  }
  
  ## Get the IPW weights, do weighted regression.
  ipw.weights <- R / ipw.scores + (1 - R) / (1 - ipw.scores)
  if (is.matrix(X)) {
    ipw.est <- summary(lm(  Y[R == 1] ~ X[R == 1, ], weights = ipw.weights[R == 1]  ))
  } else {
    ipw.est <- summary(lm(  Y[R == 1] ~ X[R == 1], weights = ipw.weights[R == 1]  ))
  }
  return(list("est" = ipw.est$coefficients, "w.min" = min(ipw.weights), "w.max" = max(ipw.weights)))
  
}

## Inverse probability weighting - logistic case.
ipw.logistic <- function (X, Y, R, Z = NA, drop.x = FALSE) {
  
  ## Regress participation on covariates, compute scores.
  if (all(is.na(Z))) {
    ipw.par <- unname(glm(R ~ X, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X) %*% ipw.par))
  } else if (drop.x == TRUE) {
    ipw.par <- unname(glm(R ~ Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, Z) %*% ipw.par))
  } else {
    ipw.par <- unname(glm(R ~ X + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X, Z) %*% ipw.par))
  }
  
  ## Get the IPW weights, do weighted regression.
  ipw.weights <- R / ipw.scores + (1 - R) / (1 - ipw.scores)
  ipw.est <- summary(glm(  Y[R == 1] ~ X[R == 1], family = binomial, weights = ipw.weights[R == 1]  ))
  return(list("est" = ipw.est$coefficients, "w.min" = min(ipw.weights), "w.max" = max(ipw.weights)))
  
}

## Inverse probability weighting - logistic case.
ipw.poisson <- function (X, Y, R, Z = NA, drop.x = FALSE) {
  
  ## Regress participation on covariates, compute scores.
  if (all(is.na(Z))) {
    ipw.par <- unname(glm(R ~ X, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X) %*% ipw.par))
  } else if (drop.x == TRUE) {
    ipw.par <- unname(glm(R ~ Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, Z) %*% ipw.par))
  } else {
    ipw.par <- unname(glm(R ~ X + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, X, Z) %*% ipw.par))
  }
  
  ## Get the IPW weights, do weighted regression.
  ipw.weights <- R / ipw.scores + (1 - R) / (1 - ipw.scores)
  ipw.est <- summary(glm(  Y[R == 1] ~ X[R == 1], family = poisson, weights = ipw.weights[R == 1]  ))
  return(list("est" = ipw.est$coefficients, "w.min" = min(ipw.weights), "w.max" = max(ipw.weights)))
  
}

## Inverse probability weighting - instrumental variables analysis.
ipw.ivreg <- function (G, X, Y, R, Z = NA, selection = "G") {
  
  ## Regress participation on covariates, compute scores.
  if (selection == "G") {
    ipw.par <- unname(glm(R ~ G, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G) %*% ipw.par))
  } else if (selection == "Z") {
    ipw.par <- unname(glm(R ~ Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, Z) %*% ipw.par))
  } else if (selection == "GZ") {
    ipw.par <- unname(glm(R ~ G + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G, Z) %*% ipw.par))
  } else if (selection == "GX") {
    ipw.par <- unname(glm(R ~ G + X, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G, X) %*% ipw.par))
  } else if (selection == "GXZ") {
    ipw.par <- unname(glm(R ~ G + X + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G, X, Z) %*% ipw.par))
  } else if (selection == "GY") {
    ipw.par <- unname(glm(R ~ G + Y, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G, Y) %*% ipw.par))
  } else if (selection == "GYZ") {
    ipw.par <- unname(glm(R ~ G + Y + Z, family = binomial)$coef)
    ipw.scores <- expit(as.vector(cbind(1, G, Y, Z) %*% ipw.par))
  } else {
    stop("Wrong selection specification.")
  }
  
  ## Get the IPW weights, do weighted IV regression.
  ipw.weights <- R / ipw.scores + (1 - R) / (1 - ipw.scores)
  if (is.matrix(G)) {
    ipw.est <- summary(ivreg(  Y[R == 1] ~ X[R == 1] | G[R == 1, ], weights = ipw.weights[R == 1]  ))
  } else {
    ipw.est <- summary(ivreg(  Y[R == 1] ~ X[R == 1] | G[R == 1], weights = ipw.weights[R == 1]  ))
  }
  return(list("est" = ipw.est$coefficients, "w.min" = min(ipw.weights), "w.max" = max(ipw.weights)))
  
}

## IV-TTW - linear case - partial likelihood.
partial.lik1 <- function (alpha, X, Z, R) {
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  sum(R * Va - log(1 + exp(Va)))
}
partial.lik2 <- function (theta, X, Z, R, Y, alpha.hat, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  log.sigma2 = theta[length(theta)]
  beta <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 3) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 2 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 2 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  
  ## Compute the likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 )
  
}

## IV-TTW - linear case - full likelihood.
full.lik <- function (theta, X, Z, R, Y, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  log.sigma2 = theta[length(theta)]
  beta <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    alpha <- theta[(2 * d + 3):(3 * d + 4)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 6) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 3 * d + 6) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 3 * d + 7) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    alpha <- theta[(2 * d + 5):(3 * d + 6)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  
  ## Compute the log-likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 ) + sum(R * Va - log(1 + exp(Va)))
  
}

## IV-TTW - logistic case - full likelihood.
full.logit.lik <- function (theta, X, Z, R, Y, Z.term = FALSE, XZ.term = FALSE, stable = TRUE) {
  
  ## Get the number of covariates and the full regression parameters.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  psi <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    alpha <- theta[(2 * d + 3):(3 * d + 4)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 3 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 3 * d + 6) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    alpha <- theta[(2 * d + 5):(3 * d + 6)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute the remaining linear terms.
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  Wp <- as.vector(cbind(1, X) %*% psi)
  
  ## Compute r.probs.
  if (stable == TRUE) {
    r.probs <- (1 - expit2(Wp)) * expit2(Va) + expit2(Wp) * expit2(Va + Wh)
  } else {
    r.probs <- (1 / (1 + exp(Wp))) * (exp(Va) / (1 + exp(Va))) + (exp(Wp) / (1 + exp(Wp))) * (exp(Va + Wh) / (1 + exp(Va + Wh)))
  }
  
  ## Compute the log-likelihood.
  y.logit <- Wp + Wh - log( (exp(Va + Wh) + 1) / (1 + exp(Va)) )
  y.log <- y.logit - log(1 + exp(y.logit))
  y.nlog <- - log(1 + exp(y.logit))
  sum(R * Y * y.log) + sum(R * (1 - Y) * y.nlog) + sum(R * log(r.probs)) + sum((1 - R) * log(1 - r.probs))
  
}

## IV-TTW - Poisson case - partial likelihood.
partial.poisson.lik1 <- function (alpha, X, Z, R) {
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  sum(R * Va - log(1 + exp(Va)))
}
partial.poisson.lik2 <- function (theta, X, Z, R, Y, alpha.hat, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Get the number of covariates and the full regression parameters.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  psi <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 2) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 3) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 2 * d + 3) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 2 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  Wp <- as.vector(cbind(1, X) %*% psi)
  
  ## Compute the likelihood.
  nu.bar <- log( (1 + exp(Wh + Va)) / (1 + exp(Va)) )
  m <- exp(Wh - nu.bar + Wp)
  sum(R * Y * log(m)) - sum(R * m)
  
}

## IV-TTW - Poisson case - full likelihood.
full.poisson.lik <- function (theta, X, Z, R, Y, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Get the number of covariates and the full regression parameters.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  psi <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    alpha <- theta[(2 * d + 3):(3 * d + 4)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 3 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(3 * d + 5)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 3 * d + 6) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    alpha <- theta[(2 * d + 5):(3 * d + 6)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  Wp <- as.vector(cbind(1, X) %*% psi)
  
  ## Compute the log-likelihood.
  nu.bar <- log( (1 + exp(Wh + Va)) / (1 + exp(Va)) )
  m <- exp(Wh - nu.bar + Wp)
  sum(R * Va - log(1 + exp(Va))) + sum(R * Y * log(m)) - sum(R * m)
  
}

## IV-TTW - linear case - partial likelihood with added variables in the selection model.
partial.lik2.v2 <- function (theta, X, Z, R, Y, alpha.hat, add.sel = NULL) {
  
  ## Identify the length of the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (is.null(add.sel)) k <- 0 else if (is.vector(add.sel)) k <- 1 else k <- ncol(add.sel)
  if (length(theta) != 2 * d + k + 3) stop("You messed up again!")
  
  ## Set up the parameter vector.
  beta <- theta[1:(d + 1)]
  eta <- theta[(d + 2):(2 * d + k + 2)]
  log.sigma2 = theta[2 * d + k + 3]
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  if (is.null(add.sel)) Wh <- as.vector(cbind(1, X) %*% eta) else Wh <- as.vector(cbind(1, X, add.sel) %*% eta) 
  Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  
  ## Compute the likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 )
  
}

## IV-TTW - linear case - full likelihood with added variables in the selection model.
full.lik.v2 <- function (theta, X, Z, R, Y, add.sel = NULL) {
  
  ## Identify the length of the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (is.vector(Z)) P <- 1 else P <- ncol(Z)
  if (is.null(add.sel)) k <- 0 else if (is.vector(add.sel)) k <- 1 else k <- ncol(add.sel)
  if (length(theta) != 3 * d + P + k + 4) stop("You messed up again!")
  
  ## Set up the parameter vector.
  beta <- theta[1:(d + 1)]
  eta <- theta[(d + 2):(2 * d + k + 2)]
  alpha <- theta[(2 * d + k + 3):(3 * d + P + k + 3)]
  log.sigma2 = theta[3 * d + P + k + 4]
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  if (is.null(add.sel)) Wh <- as.vector(cbind(1, X) %*% eta) else Wh <- as.vector(cbind(1, X, add.sel) %*% eta) 
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  
  ## Compute the log-likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 ) + sum(R * Va - log(1 + exp(Va)))
  
}

save.image("IVsel_Functions.RData")

##################################################

## The following modifications are not actually 
## included in the functions we use.


## Account for no X-R effects.

## IV-TTW - linear case - partial likelihood.
partial.lik1 <- function (alpha, X, Z, R, xr = TRUE) {
  if (xr == TRUE) {
    Va <- as.vector(cbind(1, X, Z) %*% alpha)
  } else {
    Va <- as.vector(cbind(1, Z) %*% alpha)
  }
  sum(R * Va - log(1 + exp(Va)))
}
partial.lik2 <- function (theta, X, Z, R, Y, alpha.hat, xr = TRUE, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  log.sigma2 = theta[length(theta)]
  beta <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 3) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 2 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 2 * d + 4) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 2 * d + 5) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  if (xr == TRUE) {
    Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  } else {
    Va <- as.vector(cbind(1, Z) %*% alpha.hat)
  }
  
  ## Compute the likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 )
  
}

## IV-TTW - linear case - full likelihood.
full.lik <- function (theta, X, Z, R, Y, xr = TRUE, Z.term = FALSE, XZ.term = FALSE) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  log.sigma2 = theta[length(theta)]
  beta <- theta[1:(d + 1)]
  
  ## Distinguish cases depending on which terms are fitted to delta.
  if (Z.term == FALSE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 4 + as.numeric(xr)) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 2)]
    alpha <- theta[(2 * d + 3):(length(theta) - 1)]
    Wh <- as.vector(cbind(1, X) %*% eta)
  } else if (Z.term == TRUE & XZ.term == FALSE) {
    if (length(theta) != 3 * d + 5 + as.numeric(xr)) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(length(theta) - 1)]
    Wh <- as.vector(cbind(1, X, Z) %*% eta)
  } else if (Z.term == FALSE & XZ.term == TRUE) {
    if (length(theta) != 3 * d + 5 + as.numeric(xr)) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 3)]
    alpha <- theta[(2 * d + 4):(length(theta) - 1)]
    Wh <- as.vector(cbind(1, X, X * Z) %*% eta)
  } else {
    if (length(theta) != 3 * d + 6 + as.numeric(xr)) stop("You messed up again!")
    eta <- theta[(d + 2):(2 * d + 4)]
    alpha <- theta[(2 * d + 5):(length(theta) - 1)]
    Wh <- as.vector(cbind(1, X, Z, X * Z) %*% eta)
  }
  
  ## Compute linear terms.
  Wb <- as.vector(cbind(1, X) %*% beta)
  if (xr == TRUE) {
    Va <- as.vector(cbind(1, X, Z) %*% alpha)
  } else {
    Va <- as.vector(cbind(1, Z) %*% alpha)
  }
  
  ## Compute the log-likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2 ) + sum(R * Va - log(1 + exp(Va)))
  
}
