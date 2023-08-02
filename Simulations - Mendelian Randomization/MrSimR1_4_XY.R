
##########   MR SIMULATION R1-4-XY   ##########

## Apply the IVsel methods in Mendelian randomization
## with different values for the X-Y causal effect.

## For simplicity, we run the revision MR simulations
## in a one-sample framework using a single instrument
## for inference (e.g. a PRS).

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)

## Generate array ID and file name.
range <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
filename <- paste("MrSimR1_4_XY_", range, "_results.RData", sep = "")

## Set up the values of the causal effect parameter.
beta.y.values <- c(0, 0.05, 0.1, 0.15, 0.25, 0.3)   ## 0.2 is the default, so don't repeat.

## Set up the data generation.
seed <- 3521634
iter <- 10000
n <- 10000

## Set values for the main simulation parameters.
beta.x <- sqrt(2/19)
gamma.x <- 1
beta.y <- beta.y.values[range]
gamma.y <- 1
alpha.r <- 0
beta.r <- 0.5
gamma.r <- 0.5
delta.r <- 0.5

## Make sure the proportions of genetic variation 
## in the risk factor and outcome are reasonable.
gv.x <- (beta.x^2 * 1) / (beta.x^2 + gamma.x^2 + 1)
gv.y <- (beta.y^2 * beta.x^2) / (beta.y^2 * beta.x^2 + (beta.y * gamma.x + gamma.y)^2 + beta.y^2 + 1)

## Store G-X associations here.
oracle.x <- matrix(0, iter, 4); colnames(oracle.x) <- c("Intercept", "se.int", "Beta", "se.beta")
naive.x <- oracle.x; ipw1.x <- oracle.x; ipw2.x <- oracle.x; ipw3.x <- oracle.x
heckman1.x <- oracle.x; heckman2.x <- oracle.x; heckman3.x <- oracle.x
heckman4.x <- oracle.x; partial.x <- oracle.x; full.x <- oracle.x

## Store G-Y associations here.
oracle.y <- matrix(0, iter, 4); colnames(oracle.y) <- c("Intercept", "se.int", "Beta", "se.beta")
naive.y <- oracle.y; ipw1.y <- oracle.y; ipw2.y <- oracle.y; ipw3.y <- oracle.y
heckman1.y <- oracle.y; heckman2.y <- oracle.y; heckman3.y <- oracle.y
heckman4.y <- oracle.y; partial.y <- oracle.y; full.y <- oracle.y

## Store additional diagnostics here.
diagnostics <- matrix(0, iter, 23)
colnames(diagnostics) <- c("Strength", "Strength | X", "Strength | X, Y", "Prop Selected", "MD X", "OR X", "MD Y", "OR Y", 
                           "Wmin | X", "Wmax | X", "Wmin | X, Z", "Wmax | X, Z", "MD Z", "OR Z", "F-stat G", "GV% X", "GV% Y",
                           "R2 (X ~ G)", "R2 (Y ~ G)", "McFadden Z", "McFadden X", "McFadden Y", "McFadden XY")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate the risk factor, outcome and instruments.
  U <- rnorm(n, 0, 1)
  Z <- rnorm(n, 0, 1)
  G <- rnorm(n, 0, 1)
  X <- beta.x * G + gamma.x * U + rnorm(n, 0, 1)
  Y <- beta.y * X + gamma.y * U + rnorm(n, 0, 1)
  
  ## Simulate the selection coefficient.
  R.probs <- expit(alpha.r + beta.r * X + gamma.r * Z + delta.r * Y)
  R <- rbinom(n, 1, R.probs)
  
  ## G-X association estimates.
  
  ## Naive and oracle estimates.
  naive.est.x <- summary(lm(X[R == 1] ~ G[R == 1]))
  oracle.est.x <- summary(lm(X ~ G))
  
  ## Inverse probability weighting.
  ipw.est1.x <- ipw.linear(X = G, Y = X, R = R)
  ipw.est2.x <- ipw.linear(X = G, Y = X, R = R, Z = Z)
  ipw.est3.x <- ipw.linear(X = G, Y = X, R = R, Z = Z, drop.x = TRUE)
  
  ## Heckman's method.
  Xobs <- X;  Xobs[R == 0] <- NA
  heck1.x <- summary(selection(R ~ Z, Xobs ~ G, method = "2step"))
  heck2.x <- summary(selection(R ~ Z, Xobs ~ G, method = "ml"))
  heck3.x <- summary(selection(R ~ G + Z, Xobs ~ G, method = "2step"))
  heck4.x <- summary(selection(R ~ G + Z, Xobs ~ G, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1.x <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z)
  partial.opt2.x <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = X, alpha.hat = partial.opt1.x$par, hessian = TRUE)
  partial.sd.x <- sqrt(diag(solve(- partial.opt2.x$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt.x <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = X, hessian = TRUE)
  full.sd.x <- sqrt(diag(solve(- full.opt.x$hessian))[1:2])
  
  ## G-Y association estimates.
  
  ## Naive and oracle estimates.
  naive.est.y <- summary(lm(Y[R == 1] ~ G[R == 1]))
  oracle.est.y <- summary(lm(Y ~ G))
  
  ## Inverse probability weighting.
  ipw.est1.y <- ipw.linear(X = G, Y = Y, R = R)
  ipw.est2.y <- ipw.linear(X = G, Y = Y, R = R, Z = Z)
  ipw.est3.y <- ipw.linear(X = G, Y = Y, R = R, Z = Z, drop.x = TRUE)
  
  ## Heckman's method.
  Yobs <- Y;  Yobs[R == 0] <- NA
  heck1.y <- summary(selection(R ~ Z, Yobs ~ G, method = "2step"))
  heck2.y <- summary(selection(R ~ Z, Yobs ~ G, method = "ml"))
  heck3.y <- summary(selection(R ~ G + Z, Yobs ~ G, method = "2step"))
  heck4.y <- summary(selection(R ~ G + Z, Yobs ~ G, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1.y <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z)
  partial.opt2.y <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = Y, alpha.hat = partial.opt1.y$par, hessian = TRUE)
  partial.sd.y <- sqrt(diag(solve(- partial.opt2.y$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt.y <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = Y, hessian = TRUE)
  full.sd.y <- sqrt(diag(solve(- full.opt.y$hessian))[1:2])
  
  ## Store G-X results.
  oracle.x[I, ] <- c(oracle.est.x$coef[1, 1], oracle.est.x$coef[1, 2], oracle.est.x$coef[2, 1], oracle.est.x$coef[2, 2])
  naive.x[I, ] <- c(naive.est.x$coef[1, 1], naive.est.x$coef[1, 2], naive.est.x$coef[2, 1], naive.est.x$coef[2, 2])
  ipw1.x[I, ] <- c(ipw.est1.x$est[1, 1], ipw.est1.x$est[1, 2], ipw.est1.x$est[2, 1], ipw.est1.x$est[2, 2])
  ipw2.x[I, ] <- c(ipw.est2.x$est[1, 1], ipw.est2.x$est[1, 2], ipw.est2.x$est[2, 1], ipw.est2.x$est[2, 2])
  ipw3.x[I, ] <- c(ipw.est3.x$est[1, 1], ipw.est3.x$est[1, 2], ipw.est3.x$est[2, 1], ipw.est3.x$est[2, 2])
  partial.x[I, ] <- c(partial.opt2.x$par[1], partial.sd.x[1], partial.opt2.x$par[2], partial.sd.x[2])
  full.x[I, ] <- c(full.opt.x$par[1], full.sd.x[1], full.opt.x$par[2], full.sd.x[2])
  heckman1.x[I, ] <- c(heck1.x$estimate[3, 1:2], heck1.x$estimate[4, 1:2])
  heckman2.x[I, ] <- c(heck2.x$estimate[3, 1:2], heck2.x$estimate[4, 1:2])
  heckman3.x[I, ] <- c(heck3.x$estimate[4, 1:2], heck3.x$estimate[5, 1:2])
  heckman4.x[I, ] <- c(heck4.x$estimate[4, 1:2], heck4.x$estimate[5, 1:2])
  
  ## Store G-Y results.
  oracle.y[I, ] <- c(oracle.est.y$coef[1, 1], oracle.est.y$coef[1, 2], oracle.est.y$coef[2, 1], oracle.est.y$coef[2, 2])
  naive.y[I, ] <- c(naive.est.y$coef[1, 1], naive.est.y$coef[1, 2], naive.est.y$coef[2, 1], naive.est.y$coef[2, 2])
  ipw1.y[I, ] <- c(ipw.est1.y$est[1, 1], ipw.est1.y$est[1, 2], ipw.est1.y$est[2, 1], ipw.est1.y$est[2, 2])
  ipw2.y[I, ] <- c(ipw.est2.y$est[1, 1], ipw.est2.y$est[1, 2], ipw.est2.y$est[2, 1], ipw.est2.y$est[2, 2])
  ipw3.y[I, ] <- c(ipw.est3.y$est[1, 1], ipw.est3.y$est[1, 2], ipw.est3.y$est[2, 1], ipw.est3.y$est[2, 2])
  partial.y[I, ] <- c(partial.opt2.y$par[1], partial.sd.y[1], partial.opt2.y$par[2], partial.sd.y[2])
  full.y[I, ] <- c(full.opt.y$par[1], full.sd.y[1], full.opt.y$par[2], full.sd.y[2])
  heckman1.y[I, ] <- c(heck1.y$estimate[3, 1:2], heck1.y$estimate[4, 1:2])
  heckman2.y[I, ] <- c(heck2.y$estimate[3, 1:2], heck2.y$estimate[4, 1:2])
  heckman3.y[I, ] <- c(heck3.y$estimate[4, 1:2], heck3.y$estimate[5, 1:2])
  heckman4.y[I, ] <- c(heck4.y$estimate[4, 1:2], heck4.y$estimate[5, 1:2])
  
  ## Store diagnostics.
  diagnostics[I, 1] <- glm(R ~ 1, family = binomial)$deviance - glm(R ~ Z, family = binomial)$deviance
  diagnostics[I, 2] <- glm(R ~ X, family = binomial)$deviance - glm(R ~ X + Z, family = binomial)$deviance
  diagnostics[I, 3] <- glm(R ~ X + Y, family = binomial)$deviance - glm(R ~ X + Y + Z, family = binomial)$deviance
  diagnostics[I, 4] <- mean(R)
  diagnostics[I, 5] <- mean(R[X > 0]) - mean(R[X < 0])
  diagnostics[I, 6] <- (mean(R[X > 0]) * (1 - mean(R[X < 0]))) / (mean(R[X < 0]) * (1 - mean(R[X > 0])))
  diagnostics[I, 7] <- mean(R[Y > 1]) - mean(R[Y < 1])
  diagnostics[I, 8] <- (mean(R[Y > 1]) * (1 - mean(R[Y < 1]))) / (mean(R[Y < 1]) * (1 - mean(R[Y > 1])))
  diagnostics[I, 9] <- ipw.est1.x$w.min
  diagnostics[I, 10] <- ipw.est1.x$w.max
  diagnostics[I, 11] <- ipw.est2.x$w.min
  diagnostics[I, 12] <- ipw.est2.x$w.max
  diagnostics[I, 13] <- mean(R[Z > 0]) - mean(R[Z < 0])
  diagnostics[I, 14] <- (mean(R[Z > 0]) * (1 - mean(R[Z < 0]))) / (mean(R[Z < 0]) * (1 - mean(R[Z > 0])))
  diagnostics[I, 15] <- summary(lm(X ~ G))$f[1]
  diagnostics[I, 16] <- beta.x^2 * var(G) / var(X)
  diagnostics[I, 17] <- beta.x^2 * beta.y^2 * var(G) / var(Y)
  diagnostics[I, 18] <- summary(lm(X ~ G))$r.sq
  diagnostics[I, 19] <- summary(lm(Y ~ G))$r.sq
  diagnostics[I, 20] <- 1 - glm(R ~ Z, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 21] <- 1 - glm(R ~ X, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 22] <- 1 - glm(R ~ Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 23] <- 1 - glm(R ~ X + Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  
  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 100 == 0) save(beta.x, gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n, diagnostics,
                          oracle.x, naive.x, partial.x, full.x, ipw1.x, ipw2.x, ipw3.x,
                          heckman1.x, heckman2.x, heckman3.x, heckman4.x, 
                          oracle.y, naive.y, partial.y, full.y, ipw1.y, ipw2.y, ipw3.y,
                          heckman1.y, heckman2.y, heckman3.y, heckman4.y, 
                          file = filename)
  
}
