
##########   SIMULATION S20   ##########

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)

## Results will be saved here.
filename <- "RegS20_results.RData"

## Set up the data generation.
seed <- 1068865
iter <- 10000
n <- 10000

## Set values for the main simulation parameters.
alpha.y <- 1
beta.y <- 0.1
alpha.r <- -0.5
beta.r <- 0.5
gamma.r <- 0
delta.r <- 0.5

## Store results here.
oracle <- matrix(0, iter, 4); colnames(oracle) <- c("Intercept", "se.int", "Beta", "se.beta")
naive <- oracle; ipw1 <- oracle; ipw2 <- oracle; partial <- oracle; full <- oracle
heckman1 <- oracle; heckman2 <- oracle; heckman3 <- oracle; heckman4 <- oracle; ipw3 <- oracle

## Store additional diagnostics here.
diagnostics <- matrix(0, iter, 18)
colnames(diagnostics) <- c("Strength", "Strength | X", "Strength | X, Y", "Prop Selected", "MD X", "OR X", 
                           "MD Y", "OR Y", "Wmin | X", "Wmax | X", "Wmin | X, Z", "Wmax | X, Z", "MD Z", "OR Z",
                           "McFadden Z", "McFadden X", "McFadden Y", "McFadden XY")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate the covariates, outcome and instrument.
  X <- rnorm(n, 0, 1)
  Z <- rnorm(n, 0, 1)
  Y <- alpha.y + X * beta.y + rnorm(n, 0, 1)
  
  ## Simulate the selection coefficient.
  R.probs <- expit(alpha.r + beta.r * X + gamma.r * Z + delta.r * Y)
  R <- rbinom(n, 1, R.probs)
  
  ## Naive and oracle estimates.
  naive.est <- summary(lm(Y[R == 1] ~ X[R == 1]))
  oracle.est <- summary(lm(Y ~ X))
  
  ## Inverse probability weighting.
  ipw.est1 <- ipw.linear(X, Y, R)
  ipw.est2 <- ipw.linear(X, Y, R, Z)
  ipw.est3 <- ipw.linear(X, Y, R, Z, drop.x = TRUE)
  
  ## Heckman's method.
  Yobs <- Y;  Yobs[R == 0] <- NA
  heck1 <- summary(selection(R ~ Z, Yobs ~ X, method = "2step"))
  heck2 <- summary(selection(R ~ Z, Yobs ~ X, method = "ml"))
  heck3 <- summary(selection(R ~ X + Z, Yobs ~ X, method = "2step"))
  heck4 <- summary(selection(R ~ X + Z, Yobs ~ X, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z)
  partial.opt2 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, Y = Y, alpha.hat = partial.opt1$par, hessian = TRUE)
  partial.sd <- sqrt(diag(solve(- partial.opt2$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = X, R = R, Z = Z, Y = Y, hessian = TRUE)
  full.sd <- sqrt(diag(solve(- full.opt$hessian))[1:2])
  
  ## Store results.
  oracle[I, ] <- c(oracle.est$coef[1, 1], oracle.est$coef[1, 2], oracle.est$coef[2, 1], oracle.est$coef[2, 2])
  naive[I, ] <- c(naive.est$coef[1, 1], naive.est$coef[1, 2], naive.est$coef[2, 1], naive.est$coef[2, 2])
  ipw1[I, ] <- c(ipw.est1$est[1, 1], ipw.est1$est[1, 2], ipw.est1$est[2, 1], ipw.est1$est[2, 2])
  ipw2[I, ] <- c(ipw.est2$est[1, 1], ipw.est2$est[1, 2], ipw.est2$est[2, 1], ipw.est2$est[2, 2])
  partial[I, ] <- c(partial.opt2$par[1], partial.sd[1], partial.opt2$par[2], partial.sd[2])
  full[I, ] <- c(full.opt$par[1], full.sd[1], full.opt$par[2], full.sd[2])
  
  ## Store more results.
  heckman1[I, ] <- c(heck1$estimate[3, 1:2], heck1$estimate[4, 1:2])
  heckman2[I, ] <- c(heck2$estimate[3, 1:2], heck2$estimate[4, 1:2])
  heckman3[I, ] <- c(heck3$estimate[4, 1:2], heck3$estimate[5, 1:2])
  heckman4[I, ] <- c(heck4$estimate[4, 1:2], heck4$estimate[5, 1:2])
  ipw3[I, ] <- c(ipw.est3$est[1, 1], ipw.est3$est[1, 2], ipw.est3$est[2, 1], ipw.est3$est[2, 2])
  
  ## Store diagnostics.
  diagnostics[I, 1] <- glm(R ~ 1, family = binomial)$deviance - glm(R ~ Z, family = binomial)$deviance
  diagnostics[I, 2] <- glm(R ~ X, family = binomial)$deviance - glm(R ~ X + Z, family = binomial)$deviance
  diagnostics[I, 3] <- glm(R ~ X + Y, family = binomial)$deviance - glm(R ~ X + Y + Z, family = binomial)$deviance
  diagnostics[I, 4] <- mean(R)
  diagnostics[I, 5] <- mean(R[X > 0]) - mean(R[X < 0])
  diagnostics[I, 6] <- (mean(R[X > 0]) * (1 - mean(R[X < 0]))) / (mean(R[X < 0]) * (1 - mean(R[X > 0])))
  diagnostics[I, 7] <- mean(R[Y > 1]) - mean(R[Y < 1])
  diagnostics[I, 8] <- (mean(R[Y > 1]) * (1 - mean(R[Y < 1]))) / (mean(R[Y < 1]) * (1 - mean(R[Y > 1])))
  diagnostics[I, 9] <- ipw.est1$w.min
  diagnostics[I, 10] <- ipw.est1$w.max
  diagnostics[I, 11] <- ipw.est2$w.min
  diagnostics[I, 12] <- ipw.est2$w.max
  diagnostics[I, 13] <- mean(R[Z > 0]) - mean(R[Z < 0])
  diagnostics[I, 14] <- (mean(R[Z > 0]) * (1 - mean(R[Z < 0]))) / (mean(R[Z < 0]) * (1 - mean(R[Z > 0])))
  diagnostics[I, 15] <- 1 - glm(R ~ Z, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 16] <- 1 - glm(R ~ X, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 17] <- 1 - glm(R ~ Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 18] <- 1 - glm(R ~ X + Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  
  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 100 == 0) save(alpha.y, beta.y, alpha.r, beta.r, gamma.r, delta.r, n, oracle, naive, partial, full, 
                          diagnostics, ipw1, ipw2, ipw3, heckman1, heckman2, heckman3, heckman4, 
                          file = filename)
  
}
