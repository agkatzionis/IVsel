
####################   IVSEL SOURCE CODE   ####################

## --------------------------------------------------------- ##

## This code is still under development. If you spot any bugs
## or have any suggestions on how to improve the code,
## please contact apostolos.gkatzionis@bristol.ac.uk.

## --------------------------------------------------------- ##

##################################################

## This file contains the R code to implement the "Instruments  
## for selection" methods in regression analyses and Mendelian 
## randomization studies. The relevant methods were described 
## in the paper:

## "Using Instruments for Selection to Adjust for Selection Bias 
## in Mendelian Randomization" (A. Gkatzionis et al. 2022).


## We provide functions to implement IPW and the method of
## Tchetgen Tchetgen and Wirth (2017) in a regression analysis 
## for the effects of one or more covariates X on an outcome Y. 
## We do not explicitly implement Mendelian randomization, but 
## the methods included here can be used for MR if implemented
## separately for the G-X and G-Y associations (the SNPs will 
## play the role of the "covariates" in the function calls)
## or for the first stage and the second stage of a 2SLS fit.

## We do not include implementations of Heckman's sample selection 
## model, as these are already available in the following package:
#library(sampleSelection)


## We will also make use of the following R packages.
require(MCMCpack)

##################################################

##########   PRELIMINARY FUNCTIONS   ##########

## Simple functions for expit and logit, needed later in the code.

## Compute the logit function.
logit <- function (x) log(x / (1 - x))

## Compute the expit function (this version makes it stable for x > 710).
expit <- function (x) {
  res <- exp(x) / (1 + exp(x))
  res[x > 0] <- 1 / (1 + exp(- x[x > 0]))
  return (res)
}

##################################################

##########   IPW - REGRESSION   ##########

## Perform IPW for a generalized linear model with missing outcome data.
## This fits a logistic function for missingness with no interactions.
ipw.glm <- function (X, Y, Z = NA, x.in.w = TRUE, xy.family = gaussian, sel.family = binomial(link = "logit")) {
  
  ## Ensure arguments are of the right type.
  if (is.data.frame(X)) {
    stop("We are sorry that the current implementation of ipw.glm does not admit data frames as arguments. Please convert X to a matrix.")
  } else if (is.vector(X)) {
    d <- 1
  } else if (is.matrix(X)) {
    d <- ncol(X)
  } else {
    stop("X must be a vector or matrix.")
  }
  if (!is.vector(Y)) stop("Y must be a vector.")
  if (all(is.na(Z))) {
    k <- 0
  } else if (is.data.frame(Z)) {
    stop("We are sorry that the current implementation of ipw.glm does not admit data frames as arguments. Please convert Z to a matrix.")
  } else if (is.vector(Z)) {
    k <- 1
  } else if (is.matrix(Z)) {
    k <- ncol(Z)
  } else {
    stop("Z must be a vector or matrix, if not NA.")
  }
  if (!is.logical(x.in.w)) stop("Variable x.in.w must be assigned TRUE/FALSE values.")
  
  ## Ensure arguments are of the right dimension.
  if (is.vector(X)) {
    if (length(X) != length(Y)) stop("X and Y have different length.")
  } else {
    if (nrow(X) != length(Y)) stop("X and Y have different length.")
  }
  if (!(all(is.na(Z)))) {
    if (is.vector(Z)) {
      if (length(Z) != length(Y)) stop("Z and Y have different length.")
    } else if (is.matrix(Z)) {
      if (nrow(Z) != length(Y)) stop("Z and Y have different length.")
    } 
  }
  
  ## Detect which individuals have observed outcome values.
  if (sum(is.na(Y)) == 0) warning("No missing data detected for the outcome.")
  R <- !is.na(Y)
  
  ## Fit the propensity score model.
  if (all(is.na(Z))) {
    
    ## ... without any additional variables to be used in the weighting model.
    if (all(x.in.w == TRUE)) {
      ipw.par <- unname(glm(R ~ X, family = sel.family)$coef)
      ipw.scores <- expit(as.vector(cbind(1, X) %*% ipw.par))
    } else if (all(x.in.w == FALSE)) {
      stop("No variables included in the weighting model.")
    } else {
      ipw.par <- unname(glm(R ~ X[, x.in.w], family = sel.family)$coef)
      ipw.scores <- expit(as.vector(cbind(1, X[, x.in.w]) %*% ipw.par))
    }
    
  } else {
    
    ## ... with some additional variables to be used in the weighting model.
    if (all(x.in.w == TRUE)) {
      ipw.par <- unname(glm(R ~ Z + X, family = sel.family)$coef)
      ipw.scores <- expit(as.vector(cbind(1, Z, X) %*% ipw.par))
    } else if (all(x.in.w == FALSE)) {
      ipw.par <- unname(glm(R ~ Z, family = sel.family)$coef)
      ipw.scores <- expit(as.vector(cbind(1, Z) %*% ipw.par))
    } else {
      ipw.par <- unname(glm(R ~ Z + X[, x.in.w], family = sel.family)$coef)
      ipw.scores <- expit(as.vector(cbind(1, Z, X[, x.in.w]) %*% ipw.par))
    }
    
  }
  
  ## Compute IPW weights.
  ipw.weights <- R / ipw.scores + (1 - R) / (1 - ipw.scores)
  
  ## Perform weighted regression.
  if (is.vector(X)) {
    return(summary(glm(  Y[R == 1] ~ X[R == 1], weights = ipw.weights[R == 1], family = xy.family  )))
  } else {
    return(summary(glm(  Y[R == 1] ~ X[R == 1, ], weights = ipw.weights[R == 1], family = xy.family  )))
  }
  
}

## Arguments:

## - X: a vector or matrix of covariate values for the regression
##     of interest.
## - Y: a vector of outcome values for the regression of interest. 
##     Missing values should be denoted NA.
## - Z: a vector or matrix with additional variables to be included 
##     as covariates in the weighting model, but not in the
##     regression model.
## - x.in.w: a single logical value or logical vector of the same
##     length as X. If TRUE, all covariates in X are added to the
##     weighting model; if FALSE, none of the covariates is added. If
##     a logical vector, variables assigned TRUE are included in the 
##     weighting model and variables assigned FALSE are not included.
## - xy.family: Passed to "glm(Y ~ X)", the model from which to
##     estimate the X-Y association (linear regression by default).
## - sel.family: Passed to "glm(R ~ X)", the model used to compute
##     inverse probability weights (logistic regression by default).
##     Since selection is a binary variable, using families not 
##     suited for binary outcomes will return an error.

## Note: in its current form, the function only admits matrices, 
## not data frames, as arguments for X, Z. Please convert any 
## data frames to matrices before passing them on to ipw.glm.


## Examples: 

## Linear regression.

## Seed it.
set.seed(9153)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Y <- 1 + X * 0.2 + rnorm(n, 0, 1)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(lm(Y[R == 1] ~ X[R == 1]))
oracle.est <- summary(lm(Yall ~ X))

## Inverse probability weighting.
ipw.glm(X, Y, Z)
ipw.glm(X, Y, Z, x.in.w = FALSE)
ipw.glm(X, Y, Z, sel.family = binomial(link = "probit"))
ipw.glm(X, Y, Z, sel.family = binomial(link = "cloglog"))


## Logistic regression. 

## Seed it.
set.seed(9154)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Yprob <- expit(1 + X * 0.2)
Y <- rbinom(n, 1, Yprob)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1], family = binomial))
oracle.est <- summary(glm(Yall ~ X, family = binomial))

## Inverse probability weighting.
ipw.glm(X, Y, Z, xy.family = binomial(link = "logit"))
ipw.glm(X, Y, Z, x.in.w = FALSE,  xy.family = binomial(link = "logit"))
ipw.glm(X, Y, Z, xy.family = binomial(link = "probit"))
ipw.glm(X, Y, Z, x.in.w = FALSE,  xy.family = binomial(link = "probit"))
ipw.glm(X, Y, Z, xy.family = binomial(link = "probit"), sel.family = binomial(link = "probit"))
ipw.glm(X, Y, Z, x.in.w = FALSE,  xy.family = binomial(link = "probit"), sel.family = binomial(link = "probit"))


## Poisson regression.

## Seed it.
set.seed(9155)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Y.link <- exp(1 + X * 0.2)
Y <- rpois(n, Y.link)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1], family = poisson))
oracle.est <- summary(glm(Yall ~ X, family = poisson))

## Inverse probability weighting.
ipw.glm(X, Y, Z, xy.family = poisson(link = "log"))
ipw.glm(X, Y, Z, x.in.w = FALSE,  xy.family = poisson(link = "log"))
ipw.glm(X, Y, Z)
ipw.glm(X, Y, Z, x.in.w = FALSE)


## Linear regression with many X, Z.

## Seed it.
set.seed(9156)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1.2), rnorm(n, -1, 0.8))
Z <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1))
Y <- 1 + as.vector(X %*% c(0.2, 0.1, 0.1)) + rnorm(n, 0, 1)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 +  as.vector(X %*% c(0.2, 0.1, 0.1)) +  as.vector(Z %*% c(0.3, 0.2)) + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(lm(Y[R == 1] ~ X[R == 1, ]))
oracle.est <- summary(lm(Yall ~ X))

## Inverse probability weighting.
ipw.glm(X, Y, Z)
ipw.glm(X, Y, Z, x.in.w = FALSE)

##################################################

##########   HECKMAN - REGRESSION   ##########

## The implementation of Heckman's sample selection method for regression 
## analyses can be implemented using the function "selection" from package 
## "sampleSelection". With that in mind, we did not write our own 
## implementation of the function. Note that "selection" can only fit 
## Heckman's model for linear regression or binary outcomes, so 
## implementations using a Poisson model are not feasible at the moment.

##################################################

##########   TTW - LINEAR REGRESSION   ##########

## Run the TTW method for a linear regression model with missing outcome data.
ttw.linear <- function (X, Y, Z, C = NA, partial = TRUE) {
  
  ## Ensure arguments are of the right type.
  if (is.data.frame(X)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert X to a matrix.")
  } else if (is.vector(X)) {
    d <- 1
  } else if (is.matrix(X)) {
    d <- ncol(X)
  } else {
    stop("X must be a vector or matrix.")
  }
  if (!is.vector(Y)) stop("Y must be a vector.")
  if (all(is.na(C))) {
    k <- 0
  } else if (is.data.frame(C)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert C to a matrix.")
  } else if (is.vector(C)) {
    k <- 1
  } else if (is.matrix(C) | is.data.frame(C)) {
    k <- ncol(C)
  } else {
    stop("C must be a vector or matrix, if not NA.")
  }
  if (!is.vector(Z)) stop("Z must be a vector.")
  if (!is.logical(partial)) stop("Argument ''partial'' must take a TRUE/FALSE value.")
  
  ## Ensure arguments are of the right dimension.
  if (is.vector(X)) {
    if (length(X) != length(Y)) stop("X and Y have different length.")
  } else {
    if (nrow(X) != length(Y)) stop("X and Y have different length.")
  }
  if (!(all(is.na(C)))) {
    if (is.vector(C)) {
      if (length(C) != length(Y)) stop("C and Y have different length.")
    } else if (is.matrix(C)) {
      if (nrow(C) != length(Y)) stop("C and Y have different length.")
    } 
  }
  if (length(Z) != length(Y)) stop("Z and Y have different length.")
  
  ## Detect which individuals have observed outcome values.
  if (sum(is.na(Y)) == 0) warning("No missing data detected for the outcome.")
  R <- as.numeric(!(is.na(Y)))
  
  if (partial) {
    
    ## Commence partial likelihood optimization.
    
    ## Set up starting values for the optimization.
    prm1 <- rep(0, d + 2)
    prm2 <- rep(0, 2 * d + k + 3)
    
    ## Run the optimization.
    if (all(is.na(C))) {
      opt1 <- optim(par = prm1, partial.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
      opt2 <- optim(par = prm2, partial.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = NA, Y = Y, alpha.hat = opt1$par, hessian = TRUE)
    } else {
      opt1 <- optim(par = prm1, partial.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
      opt2 <- optim(par = prm2, partial.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, alpha.hat = opt1$par, hessian = TRUE)
    }
    
    ## Invert the Hessian, compute standard errors.
    inv.hessian <- solve(- opt2$hessian)
    if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
      opt2.sd <- sqrt(diag(inv.hessian))
    } else {
      warning("The Hessian matrix obtained was not positive definite. Standard error estimates may be inaccurate.")
      lambda <- 1e-8
      while (lambda < 2) {
        inv.hessian <- solve(- opt2$hessian + lambda * length(Y) * diag(nrow(opt2$hessian)))
        if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
          opt2.sd <- sqrt(diag(inv.hessian))
          break
        } else {
          lambda <- lambda * 2
        }
      }
    }
    Estimate <- opt2$par[1:(d + k + 1)]
    StdError <- opt2.sd[1:(d + k + 1)]
    
  } else {
    
    ## Commence full likelihood optimization.
    
    ## Set up starting values for the optimization.
    prm <- rep(0, 3 * d + k + 5)
    
    ## Run the optimization.
    if (all(is.na(C))) {
      opt2 <- optim(par = prm, full.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = NA, Y = Y, hessian = TRUE)
    } else {
      opt2 <- optim(par = prm, full.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, hessian = TRUE)
    }
    
    ## Invert the Hessian, compute standard errors.
    inv.hessian <- solve(- opt2$hessian)
    if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
      opt2.sd <- sqrt(diag(inv.hessian))
    } else {
      warning("The Hessian matrix obtained was not positive definite. Standard error estimates may be inaccurate.")
      lambda <- 1e-8
      while (lambda < 2) {
        inv.hessian <- solve(- opt2$hessian + lambda * length(Y) * diag(nrow(opt2$hessian)))
        if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
          opt2.sd <- sqrt(diag(inv.hessian))
          break
        } else {
          lambda <- lambda * 2
        }
      }
    }
    Estimate <- opt2$par[1:(d + k + 1)]
    StdError <- opt2.sd[1:(d + k + 1)]
    
  }
  
  ## Return results.
  res <- cbind(Estimate, StdError)
  res.nam <- "Intercept"
  if (is.vector(X)) res.nam <- c(res.nam, "X") else if (!(is.null(colnames(X)))) res.nam <- c(res.nam, colnames(X)) else res.nam <- c(res.nam, paste("X", 1:d, sep = ""))
  if (!(all(is.na(C)))) {
    if (is.vector(C)) res.nam <- c(res.nam, "C") else if (!(is.null(colnames(C)))) res.nam <- c(res.nam, colnames(C)) else res.nam <- c(res.nam, paste("C", 1:k, sep = ""))
  }
  rownames(res) <- res.nam
  return(res)
  
}

## Arguments:

## - X: a vector or matrix of covariate values for the regression
##     of interest.
## - Y: a vector of outcome values for the regression of interest. 
##     Missing values should be denoted NA.
## - Z: a vector of values for the "instrument for selection".
## - C: a vector or matrix with additional covariates to be included
##     in the X-Y regression (but not in the selection model).
##     Can be left empty (default).
## - partial: Logical - should the algorithm implement partial 
##     likelihood optimization or full MLE? Defaults to "partial"
##     which is numerically more stable.

## Note: in its current form, the function only admits matrices, 
## not data frames, as arguments for X, C. Please convert any 
## data frames to matrices before passing them on to ttw.linear.



## Auxiliary function to compute the partial likelihood for selection.
partial.lik1 <- function (alpha, X, Z, R) {
  
  ## Compute the (partial) likelihood for selection.
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  sum(R * Va - log(1 + exp(Va)))
  
}

## Auxiliary function to compute the partial likelihood for inference.
partial.lik2 <- function (theta, X, Z, R, Y, C, alpha.hat) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (all(is.na(C))) k <- 0 else if (is.vector(C)) k <- 1 else k <- ncol(C)
  beta <- theta[1:(d + k + 1)]
  eta <- theta[(d + k + 2):(2 * d + k + 2)]
  log.sigma2 = theta[2 * d + k + 3]
  
  ## Compute linear terms.
  Wh <- as.vector(cbind(1, X) %*% eta)
  if (all(is.na(C))) Wb <- as.vector(cbind(1, X) %*% beta) else Wb <- as.vector(cbind(1, X, C) %*% beta)
  Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  
  ## Compute the likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2, na.rm = TRUE )
  
}

## Auxiliary function to compute the full likelihood.
full.lik <- function (theta, X, Z, R, Y, C) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (all(is.na(C))) k <- 0 else if (is.vector(C)) k <- 1 else k <- ncol(C)
  beta <- theta[1:(d + k + 1)]
  eta <- theta[(d + k + 2):(2 * d + k + 2)]
  alpha <- theta[(2 * d + k + 3):(3 * d + k + 4)]
  log.sigma2 = theta[3 * d + k + 5]
  
  ## Compute linear terms.
  Wh <- as.vector(cbind(1, X) %*% eta)
  if (all(is.na(C))) Wb <- as.vector(cbind(1, X) %*% beta) else Wb <- as.vector(cbind(1, X, C) %*% beta)
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  
  ## Compute the likelihood.
  - sum(R == 1) / 2 * log.sigma2 - 1 / (2 * exp(log.sigma2)) * sum( R * ( Y - Wh / (1 + exp(Va)) - Wb )^2, na.rm = TRUE ) + sum(R * Va - log(1 + exp(Va)))
  
}

## Examples:

## Simple linear regression.

## Seed it.
set.seed(9253)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Y <- 1 + X * 0.2 + rnorm(n, 0, 1)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(lm(Y[R == 1] ~ X[R == 1]))
oracle.est <- summary(lm(Yall ~ X))

## Inverse probability weighting.
ipw1 <- ipw.glm(X, Y, Z = NA)
ipw2 <- ipw.glm(X, Y, Z = Z, x.in.w = FALSE)

## TTW method.
ttw.fit1 <- ttw.linear(X = X, Y = Y, Z = Z, C = NA, partial = TRUE)
ttw.fit2 <- ttw.linear(X = X, Y = Y, Z = Z, C = NA, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2


## Linear regression with one X, many C.

## Seed it.
set.seed(9254)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y <- 1 + 0.2 * X + as.vector(C %*% c(0.2, -0.2, 0.1)) + rnorm(n, 0, 1)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(lm(Y[R == 1] ~ X[R == 1] + C[R == 1, ]))
oracle.est <- summary(lm(Yall ~ X + C))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, FALSE, FALSE, FALSE))
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE)

## TTW method.
ttw.fit1 <- ttw.linear(X = X, Y = Y, Z = Z, C = C, partial = TRUE)
ttw.fit2 <- ttw.linear(X = X, Y = Y, Z = Z, C = C, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2


## Linear regression with many X, many C.

## Seed it.
set.seed(9255)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- cbind(rnorm(n, 0, 1), rnorm(n, -2, 1.2), rnorm(n, 2, 1.1))
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y <- 1 + as.vector(X %*% c(0.2, -0.1, -0.1)) + as.vector(C %*% c(0.2, -0.2, 0.1)) + rnorm(n, 0, 1)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + as.vector(X %*% c(0.2, 0.2, 0.1)) +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(lm(Y[R == 1] ~ X[R == 1, ] + C[R == 1, ]))
oracle.est <- summary(lm(Yall ~ X + C))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE)

## TTW method.
ttw.fit1 <- ttw.linear(X = X, Y = Y, Z = Z, C = C, partial = TRUE)
ttw.fit2 <- ttw.linear(X = X, Y = Y, Z = Z, C = C, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2


##################################################

##########   TTW - LOGISTIC REGRESSION   ##########

## Run the TTW method for a logistic regression model with missing outcome data.
## This only admits full-likelihood optimization (no partial likelihood).
ttw.logistic <- function (X, Y, Z, C = NA) {
  
  ## Ensure arguments are of the right type.
  if (is.data.frame(X)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert X to a matrix.")
  } else if (is.vector(X)) {
    d <- 1
  } else if (is.matrix(X)) {
    d <- ncol(X)
  } else {
    stop("X must be a vector or matrix.")
  }
  if (!is.vector(Y)) stop("Y must be a vector.")
  if (all(is.na(C))) {
    k <- 0
  } else if (is.data.frame(C)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert C to a matrix.")
  } else if (is.vector(C)) {
    k <- 1
  } else if (is.matrix(C) | is.data.frame(C)) {
    k <- ncol(C)
  } else {
    stop("C must be a vector or matrix, if not NA.")
  }
  if (!is.vector(Z)) stop("Z must be a vector.")

  ## Ensure arguments are of the right dimension.
  if (is.vector(X)) {
    if (length(X) != length(Y)) stop("X and Y have different length.")
  } else {
    if (nrow(X) != length(Y)) stop("X and Y have different length.")
  }
  if (!(all(is.na(C)))) {
    if (is.vector(C)) {
      if (length(C) != length(Y)) stop("C and Y have different length.")
    } else if (is.matrix(C)) {
      if (nrow(C) != length(Y)) stop("C and Y have different length.")
    } 
  }
  if (length(Z) != length(Y)) stop("Z and Y have different length.")
  
  ## Ensure the outcome is indeed binary.
  if (length(unique(Y[which(!is.na(Y))])) > 2) stop("Y is not a binary variable.")
  
  ## Detect which individuals have observed outcome values.
  if (sum(is.na(Y)) == 0) warning("No missing data detected for the outcome.")
  R <- as.numeric(!(is.na(Y)))
  
  ## Commence full likelihood optimization.
  
  ## Set up starting values for the optimization.
  prm <- rep(0, 3 * d + k + 4)
  
  ## Run the optimization.
  opt2 <- optim(par = prm, full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, hessian = TRUE)
  #if (all(is.na(C))) {
  #  opt2 <- optim(par = prm, full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = NA, Y = Y, hessian = TRUE)
  #} else {
  #  opt2 <- optim(par = prm, full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, hessian = TRUE)
  #}
  
  ## Invert the Hessian, compute standard errors.
  inv.hessian <- solve(- opt2$hessian)
  if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
    opt2.sd <- sqrt(diag(inv.hessian))
  } else {
    warning("The Hessian matrix obtained was not positive definite. Standard error estimates may be inaccurate.")
    lambda <- 1e-8
    while (lambda < 2) {
      inv.hessian <- solve(- opt2$hessian + lambda * length(Y) * diag(nrow(opt2$hessian)))
      if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
        opt2.sd <- sqrt(diag(inv.hessian))
        break
      } else {
        lambda <- lambda * 2
      }
    }
  }
  Estimate <- opt2$par[1:(d + k + 1)]
  StdError <- opt2.sd[1:(d + k + 1)]
  
  ## Return results.
  res <- cbind(Estimate, StdError)
  res.nam <- "Intercept"
  if (is.vector(X)) res.nam <- c(res.nam, "X") else if (!(is.null(colnames(X)))) res.nam <- c(res.nam, colnames(X)) else res.nam <- c(res.nam, paste("X", 1:d, sep = ""))
  if (!(all(is.na(C)))) {
    if (is.vector(C)) res.nam <- c(res.nam, "C") else if (!(is.null(colnames(C)))) res.nam <- c(res.nam, colnames(C)) else res.nam <- c(res.nam, paste("C", 1:k, sep = ""))
  }
  rownames(res) <- res.nam
  return(res)
  
}

## Arguments:

## - X: a vector or matrix of covariate values for the regression
##     of interest.
## - Y: a vector of outcome values for the regression of interest. 
##     Missing values should be denoted NA.
## - Z: a vector of values for the "instrument for selection".
## - C: a vector or matrix with additional covariates to be included
##     in the X-Y regression (but not in the selection model).
##     Can be left empty (default).

## Note: in its current form, the function only admits matrices, 
## not data frames, as arguments for X, C. Please convert any 
## data frames to matrices before passing them on to ttw.logistic.



## Auxiliary function that computes the (full) likelihood.
full.logit.lik <- function (theta, X, Z, R, Y, C) {
  
  ## Get the number of covariates, split the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (all(is.na(C))) k <- 0 else if (is.vector(C)) k <- 1 else k <- ncol(C)
  psi <- theta[1:(d + k + 1)]
  eta <- theta[(d + k + 2):(2 * d + k + 2)]
  alpha <- theta[(2 * d + k + 3):(3 * d + k + 4)]
  
  ## Compute basic quantities.
  Wh <- as.vector(cbind(1, X) %*% eta)
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  if (all(is.na(C))) Wp <- as.vector(cbind(1, X) %*% psi) else Wp <- as.vector(cbind(1, X, C) %*% psi)
  
  ## Estimate missingness probabilities.
  r.probs <- (1 - expit(Wp)) * expit(Va) + expit(Wp) * expit(Va + Wh)
  
  ## Compute the log-likelihood.
  y.logit <- Wp + Wh - log( (exp(Va + Wh) + 1) / (1 + exp(Va)) )
  y.log <- y.logit - log(1 + exp(y.logit))
  y.nlog <- - log(1 + exp(y.logit))
  sum(R * Y * y.log, na.rm = TRUE) + sum(R * (1 - Y) * y.nlog, na.rm = TRUE) + sum(R * log(r.probs)) + sum((1 - R) * log(1 - r.probs))
  
}


## Examples:

## Simple logistic regression.

## Seed it.
set.seed(9353)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Y.prob <- expit(1 + X * 0.2)
Y <- rbinom(n, 1, Y.prob)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1], family = binomial))
oracle.est <- summary(glm(Yall ~ X, family = binomial))

## Inverse probability weighting.
ipw1 <- ipw.glm(X, Y, Z = NA, xy.family = binomial)
ipw2 <- ipw.glm(X, Y, Z = Z, x.in.w = FALSE, xy.family = binomial)

## TTW method.
ttw.fit1 <- ttw.logistic(X = X, Y = Y, Z = Z, C = NA)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1


## Linear regression with one X, many C.

## Seed it.
set.seed(9354)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y.prob <- expit(1 + 0.2 * X + as.vector(C %*% c(0.2, -0.2, 0.1)))
Y <- rbinom(n, 1, Y.prob)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1] + C[R == 1, ], family = binomial))
oracle.est <- summary(glm(Yall ~ X + C, family = binomial))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, FALSE, FALSE, FALSE), xy.family = binomial)
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE, xy.family = binomial)

## TTW method.
ttw.fit1 <- ttw.logistic(X = X, Y = Y, Z = Z, C = C)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1


## Linear regression with many X, many C.

## Seed it.
set.seed(9355)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- cbind(rnorm(n, 0, 1), rnorm(n, -2, 1.2), rnorm(n, 2, 1.1))
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y.prob <- expit(1 + as.vector(X %*% c(0.2, -0.1, -0.1)) + as.vector(C %*% c(0.2, -0.2, 0.1)))
Y <- rbinom(n, 1, Y.prob)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + as.vector(X %*% c(0.2, 0.2, 0.1)) +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1, ] + C[R == 1, ], family = binomial))
oracle.est <- summary(glm(Yall ~ X + C, family = binomial))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), xy.family = binomial)
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE, xy.family = binomial)

## TTW method.
ttw.fit1 <- ttw.logistic(X = X, Y = Y, Z = Z, C = C)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1

##################################################

##########   TTW - POISSON REGRESSION   ##########

## Run the TTW method for a Poisson regression model with missing outcome data.
ttw.poisson <- function (X, Y, Z, C = NA, partial = TRUE) {
  
  ## Ensure arguments are of the right type.
  if (is.data.frame(X)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert X to a matrix.")
  } else if (is.vector(X)) {
    d <- 1
  } else if (is.matrix(X)) {
    d <- ncol(X)
  } else {
    stop("X must be a vector or matrix.")
  }
  if (!is.vector(Y)) stop("Y must be a vector.")
  if (all(is.na(C))) {
    k <- 0
  } else if (is.data.frame(C)) {
    stop("We are sorry that the current implementation of ttw.linear does not admit data frames as arguments. Please convert C to a matrix.")
  } else if (is.vector(C)) {
    k <- 1
  } else if (is.matrix(C) | is.data.frame(C)) {
    k <- ncol(C)
  } else {
    stop("C must be a vector or matrix, if not NA.")
  }
  if (!is.vector(Z)) stop("Z must be a vector.")
  if (!is.logical(partial)) stop("Argument ''partial'' must take a TRUE/FALSE value.")
  
  ## Ensure arguments are of the right dimension.
  if (is.vector(X)) {
    if (length(X) != length(Y)) stop("X and Y have different length.")
  } else {
    if (nrow(X) != length(Y)) stop("X and Y have different length.")
  }
  if (!(all(is.na(C)))) {
    if (is.vector(C)) {
      if (length(C) != length(Y)) stop("C and Y have different length.")
    } else if (is.matrix(C)) {
      if (nrow(C) != length(Y)) stop("C and Y have different length.")
    } 
  }
  if (length(Z) != length(Y)) stop("Z and Y have different length.")
  
  ## Check if the outcome contains integer values.
  if (!(all(Y[which(!(is.na(Y)))] - floor(Y[which(!(is.na(Y)))]) == 0))) warning("Non-integer values detected for the outcome.")
  
  ## Detect which individuals have observed outcome values.
  if (sum(is.na(Y)) == 0) warning("No missing data detected for the outcome.")
  R <- as.numeric(!(is.na(Y)))
  
  if (partial) {
    
    ## Commence partial likelihood optimization.
    
    ## Set up starting values for the optimization.
    prm1 <- rep(0, d + 2)
    prm2 <- rep(0, 2 * d + k + 2)
    
    ## Run the optimization.
    #if (all(is.na(C))) {
    #  opt1 <- optim(par = prm1, partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
    #  opt2 <- optim(par = prm2, partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = NA, Y = Y, alpha.hat = opt1$par, hessian = TRUE)
    #} else {
    #  opt1 <- optim(par = prm1, partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
    #  opt2 <- optim(par = prm2, partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, alpha.hat = opt1$par, hessian = TRUE)
    #}
    opt1 <- optim(par = prm1, partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
    opt2 <- optim(par = prm2, partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, alpha.hat = opt1$par, hessian = TRUE)
    
    ## Invert the Hessian, compute standard errors.
    inv.hessian <- solve(- opt2$hessian)
    if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
      opt2.sd <- sqrt(diag(inv.hessian))
    } else {
      warning("The Hessian matrix obtained was not positive definite. Standard error estimates may be inaccurate.")
      lambda <- 1e-8
      while (lambda < 2) {
        inv.hessian <- solve(- opt2$hessian + lambda * length(Y) * diag(nrow(opt2$hessian)))
        if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
          opt2.sd <- sqrt(diag(inv.hessian))
          break
        } else {
          lambda <- lambda * 2
        }
      }
    }
    Estimate <- opt2$par[1:(d + k + 1)]
    StdError <- opt2.sd[1:(d + k + 1)]
    
  } else {
    
    ## Commence full likelihood optimization.
    
    ## Set up starting values for the optimization.
    prm <- rep(0, 3 * d + k + 4)
    
    ## Run the optimization.
    #if (all(is.na(C))) {
    #  opt2 <- optim(par = prm, full.poisson.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = NA, Y = Y, hessian = TRUE)
    #} else {
    #  opt2 <- optim(par = prm, full.poisson.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, hessian = TRUE)
    #}
    opt2 <- optim(par = prm, full.poisson.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, C = C, Y = Y, hessian = TRUE)
    
    ## Invert the Hessian, compute standard errors.
    inv.hessian <- solve(- opt2$hessian)
    if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
      opt2.sd <- sqrt(diag(inv.hessian))
    } else {
      warning("The Hessian matrix obtained was not positive definite. Standard error estimates may be inaccurate.")
      lambda <- 1e-8
      while (lambda < 2) {
        inv.hessian <- solve(- opt2$hessian + lambda * length(Y) * diag(nrow(opt2$hessian)))
        if (all(diag(inv.hessian)[1:(d + k + 1)] > 0)) {
          opt2.sd <- sqrt(diag(inv.hessian))
          break
        } else {
          lambda <- lambda * 2
        }
      }
    }
    Estimate <- opt2$par[1:(d + k + 1)]
    StdError <- opt2.sd[1:(d + k + 1)]
    
  }
  
  ## Return results.
  res <- cbind(Estimate, StdError)
  res.nam <- "Intercept"
  if (is.vector(X)) res.nam <- c(res.nam, "X") else if (!(is.null(colnames(X)))) res.nam <- c(res.nam, colnames(X)) else res.nam <- c(res.nam, paste("X", 1:d, sep = ""))
  if (!(all(is.na(C)))) {
    if (is.vector(C)) res.nam <- c(res.nam, "C") else if (!(is.null(colnames(C)))) res.nam <- c(res.nam, colnames(C)) else res.nam <- c(res.nam, paste("C", 1:k, sep = ""))
  }
  rownames(res) <- res.nam
  return(res)
  
}

## Arguments:

## - X: a vector or matrix of covariate values for the regression
##     of interest.
## - Y: a vector of outcome values for the regression of interest. 
##     Missing values should be denoted NA.
## - Z: a vector of values for the "instrument for selection".
## - C: a vector or matrix with additional covariates to be included
##     in the X-Y regression (but not in the selection model).
##     Can be left empty (default).
## - partial: Logical - should the algorithm implement partial 
##     likelihood optimization or full MLE? Defaults to "partial"
##     which is numerically more stable.

## Note: in its current form, the function only admits matrices, 
## not data frames, as arguments for X, C. Please convert any 
## data frames to matrices before passing them on to ttw.poisson.



## Auxiliary function to compute the partial likelihood for selection.
partial.poisson.lik1 <- function (alpha, X, Z, R) {
  
  ## Compute the (partial) likelihood for selection.
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  sum(R * Va - log(1 + exp(Va)))
  
}

## Auxiliary function to compute the partial likelihood for inference.
partial.poisson.lik2 <- function (theta, X, Z, R, Y, C, alpha.hat) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (all(is.na(C))) k <- 0 else if (is.vector(C)) k <- 1 else k <- ncol(C)
  psi <- theta[1:(d + k + 1)]
  eta <- theta[(d + k + 2):(2 * d + k + 2)]

  ## Compute linear terms.
  Wh <- as.vector(cbind(1, X) %*% eta)
  if (all(is.na(C))) Wp <- as.vector(cbind(1, X) %*% psi) else Wp <- as.vector(cbind(1, X, C) %*% psi)
  Va <- as.vector(cbind(1, X, Z) %*% alpha.hat)
  
  ## Compute the likelihood.
  nu.bar <- log( (1 + exp(Wh + Va)) / (1 + exp(Va)) )
  m <- exp(Wh - nu.bar + Wp)
  sum(R * Y * log(m), na.rm = TRUE) - sum(R * m)
  
}

## Auxiliary function to compute the full likelihood.
full.poisson.lik <- function (theta, X, Z, R, Y, C) {
  
  ## Set up the parameter vector.
  if (is.vector(X)) d <- 1 else d <- ncol(X)
  if (all(is.na(C))) k <- 0 else if (is.vector(C)) k <- 1 else k <- ncol(C)
  psi <- theta[1:(d + k + 1)]
  eta <- theta[(d + k + 2):(2 * d + k + 2)]
  alpha <- theta[(2 * d + k + 3):(3 * d + k + 4)]

  ## Compute linear terms.
  Wh <- as.vector(cbind(1, X) %*% eta)
  if (all(is.na(C))) Wp <- as.vector(cbind(1, X) %*% psi) else Wp <- as.vector(cbind(1, X, C) %*% psi)
  Va <- as.vector(cbind(1, X, Z) %*% alpha)
  
  ## Compute the likelihood.
  nu.bar <- log( (1 + exp(Wh + Va)) / (1 + exp(Va)) )
  m <- exp(Wh - nu.bar + Wp)
  sum(R * Va - log(1 + exp(Va))) + sum(R * Y * log(m), na.rm = TRUE) - sum(R * m)
  
}


## Examples:

## Simple Poisson regression.

## Seed it.
set.seed(9453)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
Y.link <- exp(1 + X * 0.2)
Y <- rpois(n, Y.link)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X + 0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1], family = poisson))
oracle.est <- summary(glm(Yall ~ X, family = poisson))

## Inverse probability weighting.
ipw1 <- ipw.glm(X, Y, Z = NA, xy.family = poisson)
ipw2 <- ipw.glm(X, Y, Z = Z, x.in.w = FALSE, xy.family = poisson)

## TTW method.
ttw.fit1 <- ttw.poisson(X = X, Y = Y, Z = Z, C = NA, partial = TRUE)
ttw.fit2 <- ttw.poisson(X = X, Y = Y, Z = Z, C = NA, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2


## Linear regression with one X, many C.

## Seed it.
set.seed(9454)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- rnorm(n, 0, 1)
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y.link <- exp(1 + 0.2 * X + as.vector(C %*% c(0.2, -0.2, 0.1)))
Y <- rpois(n, Y.link)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + 0.5 * X +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1] + C[R == 1, ], family = poisson))
oracle.est <- summary(glm(Yall ~ X + C, family = poisson))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, FALSE, FALSE, FALSE), xy.family = poisson)
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE, xy.family = poisson)

## TTW method.
ttw.fit1 <- ttw.poisson(X = X, Y = Y, Z = Z, C = C, partial = TRUE)
ttw.fit2 <- ttw.poisson(X = X, Y = Y, Z = Z, C = C, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2


## Linear regression with many X, many C.

## Seed it.
set.seed(9455)

## Simulate the covariates, outcome and instrument.
n <- 10000
X <- cbind(rnorm(n, 0, 1), rnorm(n, -2, 1.2), rnorm(n, 2, 1.1))
Z <- rnorm(n, 0, 1)
C <- cbind(rnorm(n, 0, 1), rnorm(n, 2, 1), rnorm(n, -1, 0.8))
Y.link <- exp(1 + as.vector(X %*% c(0.2, -0.1, -0.1)) + as.vector(C %*% c(0.2, -0.2, 0.1)))
Y <- rpois(n, Y.link)

## Simulate the selection coefficient.
R.probs <- expit(-0.5 + as.vector(X %*% c(0.2, 0.2, 0.1)) +  0.4 * Z + 0.5 * Y)
R <- rbinom(n, 1, R.probs)
Yall <- Y
Y[R == 0] <- NA

## Naive and oracle estimates.
naive.est <- summary(glm(Y[R == 1] ~ X[R == 1, ] + C[R == 1, ], family = poisson))
oracle.est <- summary(glm(Yall ~ X + C, family = poisson))

## Inverse probability weighting.
ipw1 <- ipw.glm(cbind(X, C), Y, Z = NA, x.in.w = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), xy.family = poisson)
ipw2 <- ipw.glm(cbind(X, C), Y, Z = Z, x.in.w = FALSE, xy.family = poisson)

## TTW method.
ttw.fit1 <- ttw.poisson(X = X, Y = Y, Z = Z, C = C, partial = TRUE)
ttw.fit2 <- ttw.poisson(X = X, Y = Y, Z = Z, C = C, partial = FALSE)

## Results.
naive.est$coefficients
oracle.est$coefficients
ipw1$coefficients
ipw2$coefficients
ttw.fit1
ttw.fit2

##################################################

##########   BOUNDS FOR SELECTION BIAS   ##########

## Finally, we provide some R code to compute the selection bias
## bounds developed in Marden et al. (2018).

## The function computes nonparametric, Robins-Manski, Bayesian 
## and Asymptotic bounds, as described in Marden et al. (2018),
## and conducts a chi-squared test of association between Z and Y.
sel.bias.bounds <- function (data, alpha = 0.95, n.draws = 1000, seed = NULL) {
  
  ## Diagnostics.
  if (is.matrix(data)) {
    if (!(all(c("Total", "Missing", "Healthy", "Diseased", "Observed") %in% colnames(data)))) stop("Please specify appropriate column names for \"data\".")
  } else if (is.data.frame(data)) {
    if (!(all(c("Total", "Missing", "Healthy", "Diseased", "Observed") %in% names(data)))) stop("Please specify appropriate column names for \"data\".")
  } else {
    stop("Argument \"data\" must be a matrix or data frame.")
  }
  if (!is.numeric(alpha) | length(alpha) != 1) stop("Argument \"alpha\" must be a single number between 0 and 1.")
  if (!is.numeric(n.draws) | length(n.draws) != 1) stop("Argument \"n.draws\" must be a single number between 0 and 1.")
  if (!is.numeric(seed) | length(seed) != 1) stop("Argument \"seed\" must be a single number between 0 and 1.")
  if (!is.integer(n.draws)) n.draws <- as.integer(n.draws)
  
  ## Make sure the numbers add up.
  if (!(all(data$Diseased + data$Healthy == data$Observed))) stop("The \"Diseased\" and \"Healthy\" counts do not sum to \"Observed\".")
  if (!(all(data$Observed + data$Missing == data$Total))) stop("The \"Observed\" and \"Missing\" counts do not sum to \"Total\".")
  
  ## Set the seed, if provided.
  if (!(is.null(seed))) set.seed(seed)
  
  ## Number of IV categories.
  d <- nrow(data)
  
  ## Nonparametric bounds. 
  p1 <- sum(data$Missing) / (sum(data$Missing) + sum(data$Observed))
  p2 <- sum(data$Observed) / (sum(data$Missing) + sum(data$Observed))
  p3 <- sum(data$Diseased) / (sum(data$Diseased) + sum(data$Healthy))
  NP.bounds <- c(p3 * p2, p3 * p2 + p1)
  
  ## Robins-Manski bounds.
  p1z <- data$Missing / (data$Missing + data$Observed)
  p2z <- data$Observed / (data$Missing + data$Observed)
  p3z <- data$Diseased / (data$Diseased + data$Healthy)
  RM.bounds <- c(max(p3z * p2z), min(p3z * p2z + p1z))
  
  ## Chi-squared test.
  chisq.expected1 <- sum(data$Diseased) / sum(data$Observed) * data$Observed
  chisq.expected0 <- sum(data$Healthy) / sum(data$Observed) * data$Observed
  chisq.stat <- sum((data$Diseased - chisq.expected1)^2 / chisq.expected1) + sum((data$Healthy - chisq.expected0)^2 / chisq.expected0)
  chisq.pval <- 1 - pchisq(chisq.stat, df = 29)
  
  ## Bayesian credible sets.
  z.draws <- array(0, dim = c(n.draws, 3, d))
  for (i in 1:d) {
    z.draws[, , i] <- rdirichlet(n.draws, c(1 + data$Diseased[i], 1 + data$Healthy[i], 1 + data$Missing[i]))
  }
  z.bounds <- matrix(0, n.draws, 2)
  for (i in 1:n.draws) {
    z.bounds[i, 1] <- max(z.draws[i, 1, ])
    z.bounds[i, 2] <- min(z.draws[i, 1, ] + z.draws[i, 3, ])
  }
  prop.valid <- sum(NP.bounds[1] < z.bounds[, 1] & z.bounds[, 1] < z.bounds[, 2] & z.bounds[, 2] < NP.bounds[2]) / n.draws
  which.valid <- which(NP.bounds[1] < z.bounds[, 1] & z.bounds[, 1] < z.bounds[, 2] & z.bounds[, 2] < NP.bounds[2])
  Bayes.bounds <- c(median(z.bounds[which.valid, 1]), median(z.bounds[which.valid, 2]))
  Bayes.lower.ci <- c(sort(z.bounds[which.valid, 1])[round(prop.valid * n.draws * (1 - alpha) / 2)], sort(z.bounds[which.valid, 1])[round(prop.valid * n.draws * (1 - (1 - alpha) / 2))])
  Bayes.upper.ci <- c(sort(z.bounds[which.valid, 2])[round(prop.valid * n.draws * (1 - alpha) / 2)], sort(z.bounds[which.valid, 2])[round(prop.valid * n.draws * (1 - (1 - alpha) / 2))])
  
  
  ## Asymptotic confidence intervals.
  p.diseased <- data$Diseased / data$Total
  p.not.healthy <- (data$Diseased + data$Missing) / data$Total
  s.diseased <- sqrt(p.diseased * (1 - p.diseased) / data$Total)
  s.not.healthy <- sqrt(p.not.healthy * (1 - p.not.healthy) / data$Total)
  grid <- 10:990 * 0.001
  phi.lower <- matrix(0, 981, d)
  phi.upper <- matrix(0, 981, d)
  for (i in 1:d) {
    phi.lower[, i] <- pnorm((grid - p.diseased[i]) / s.diseased[i])
    phi.upper[, i] <- 1 - pnorm((grid - p.not.healthy[i]) / s.not.healthy[i])
  }
  phi.prod.lower <- apply(phi.lower, 1, prod)
  phi.prod.upper <- apply(phi.upper, 1, prod)
  Asymptotic.lower.ci <- c(grid[which.min(abs(phi.prod.lower - (1 - alpha) / 2))], grid[which.min(abs(phi.prod.lower - (1 - (1 - alpha) / 2)))])
  Asymptotic.upper.ci <- c(grid[which.min(abs(phi.prod.upper - (1 - (1 - alpha) / 2)))], grid[which.min(abs(phi.prod.upper - (1 - alpha) / 2))])
  Asymptotic.bounds <- c(Asymptotic.lower.ci[1], Asymptotic.upper.ci[2])
  
  ## Print outputs.
  return(list("Nonparametric bounds" = NP.bounds, "Robins-Manski bounds" = RM.bounds, 
              "Bayesian lower CI" = Bayes.lower.ci, "Bayesian upper CI" = Bayes.upper.ci, 
              "Asymptotic lower CI" = Asymptotic.lower.ci, "Asymptotic upper CI" = Asymptotic.upper.ci, 
              "Chi Squared test statistic" = chisq.stat, "Chi Squared P-value" = chisq.pval))
  
}

## Arguments:

## - data: a matrix or data frame whose rows represent different values of
##     the instrumental variable z. It should contain columns with names:
##     "Total": Total number of individuals with Z = z.
##     "Observed": Number of individuals with an observed outcome, for Z = z.
##     "Missing": Number of individuals with a missing outcome, for Z = z.
##     "Diseased": Number of individuals who are observed to have the disease, for Z = z.
##     "Healthy": Number of individuals who are observed to be healthy, for Z = z.
## - alpha: significance level for the Bayesian and asymptotic bounds.
## - n.draws: the number of random draws for the Bayesian set.
## - seed: a seed to initialize the random sampling.

## Examples:




##################################################






##################################################
