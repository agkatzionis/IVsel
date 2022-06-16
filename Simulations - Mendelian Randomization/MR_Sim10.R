
##########   MR SIMULATION 10   ##########

## One Instrument - Two Samples.

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)

## Results will be saved here.
filename <- "MR10_results.RData"

## Set up the data generation.
seed <- 67208383
iter <- 10000
n1 <- 10000
n2 <- 10000

## Set values for the main simulation parameters.
beta.x <- sqrt(2 / 19)
gamma.x <- 1
beta.y <- 0   ## Null causal effect.
gamma.y <- 1
alpha.r <- 0
beta.r <- 0   ## No X-R effect.
gamma.r <- 0.5
delta.r <- 1

## Make sure the proportions of genetic variation 
## in the risk factor and outcome are reasonable.
gv.x <- (beta.x^2 * 1) / (beta.x^2 + gamma.x^2 + 1)
gv.y <- (beta.y^2 * beta.x^2) / (beta.y^2 * beta.x^2 + (beta.y * gamma.x + gamma.y)^2 + beta.y^2 + 1)

## Store G-X associations here.
gx <- matrix(0, iter, 4); colnames(gx) <- c("Intercept", "se.int", "Beta", "se.beta")

## Store G-Y associations here.
oracle <- matrix(0, iter, 4); colnames(oracle) <- c("Intercept", "se.int", "Beta", "se.beta")
naive <- oracle; ipw1 <- oracle; ipw2 <- oracle; ipw3 <- oracle
heckman1 <- oracle; heckman2 <- oracle; heckman3 <- oracle
heckman4 <- oracle; partial <- oracle; full <- oracle

## Store additional diagnostics here.
diagnostics <- matrix(0, iter, 23)
colnames(diagnostics) <- c("Strength", "Strength | X", "Strength | X, Y", "Prop Selected", "MD X", "OR X", "MD Y", "OR Y", 
                           "Wmin | X", "Wmax | X", "Wmin | X, Z", "Wmax | X, Z", "MD Z", "OR Z", "F-stat G", "GV% X", "GV% Y",
                           "R2 (X ~ G)", "R2 (Y ~ G)", "McFadden Z", "McFadden X", "McFadden Y", "McFadden XY")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate the first sample.
  U1 <- rnorm(n1, 0, 1)
  Z1 <- rnorm(n1, 0, 1)
  G1 <- rnorm(n1, 0, 1)
  X1 <- beta.x * G1 + gamma.x * U1 + rnorm(n1, 0, 1)
  
  ## Simulate the second sample.
  U2 <- rnorm(n2, 0, 1)
  Z2 <- rnorm(n2, 0, 1)
  G2 <- rnorm(n2, 0, 1)
  X2 <- beta.x * G2 + gamma.x * U2 + rnorm(n2, 0, 1)
  Y2 <- beta.y * X2 + gamma.y * U2 + rnorm(n2, 0, 1)
  
  ## Simulate the selection coefficient.
  R2.probs <- expit(alpha.r + beta.r * X2 + gamma.r * Z2 + delta.r * Y2)
  R2 <- rbinom(n2, 1, R2.probs)
  
  ## G-X association estimates.
  gx.est <- summary(lm(X1 ~ G1))
  
  ## Naive and oracle estimates.
  naive.est <- summary(lm(Y2[R2 == 1] ~ G2[R2 == 1]))
  oracle.est <- summary(lm(Y2 ~ G2))
  
  ## Inverse probability weighting.
  ipw.est1 <- ipw.linear(X = G2, Y = Y2, R = R2)
  ipw.est2 <- ipw.linear(X = G2, Y = Y2, R = R2, Z = Z2)
  ipw.est3 <- ipw.linear(X = G2, Y = Y2, R = R2, Z = Z2, drop.x = TRUE)
  
  ## Heckman's method.
  Y2obs <- Y2;  Y2obs[R2 == 0] <- NA
  heck1 <- summary(selection(R2 ~ Z2, Y2obs ~ G2, method = "2step"))
  heck2 <- summary(selection(R2 ~ Z2, Y2obs ~ G2, method = "ml"))
  heck3 <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2, method = "2step"))
  heck4 <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2)
  partial.opt2 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2, Y = Y2, alpha.hat = partial.opt1$par, hessian = TRUE)
  partial.sd <- sqrt(diag(solve(- partial.opt2$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2, Y = Y2, hessian = TRUE)
  full.sd <- sqrt(diag(solve(- full.opt$hessian))[1:2])
  
  ## Store G-X results.
  gx[I, ] <- c(gx.est$coef[1, 1:2], gx.est$coef[2, 1:2])
  
  ## Store G-Y results.
  oracle[I, ] <- c(oracle.est$coef[1, 1], oracle.est$coef[1, 2], oracle.est$coef[2, 1], oracle.est$coef[2, 2])
  naive[I, ] <- c(naive.est$coef[1, 1], naive.est$coef[1, 2], naive.est$coef[2, 1], naive.est$coef[2, 2])
  ipw1[I, ] <- c(ipw.est1$est[1, 1], ipw.est1$est[1, 2], ipw.est1$est[2, 1], ipw.est1$est[2, 2])
  ipw2[I, ] <- c(ipw.est2$est[1, 1], ipw.est2$est[1, 2], ipw.est2$est[2, 1], ipw.est2$est[2, 2])
  ipw3[I, ] <- c(ipw.est3$est[1, 1], ipw.est3$est[1, 2], ipw.est3$est[2, 1], ipw.est3$est[2, 2])
  partial[I, ] <- c(partial.opt2$par[1], partial.sd[1], partial.opt2$par[2], partial.sd[2])
  full[I, ] <- c(full.opt$par[1], full.sd[1], full.opt$par[2], full.sd[2])
  heckman1[I, ] <- c(heck1$estimate[3, 1:2], heck1$estimate[4, 1:2])
  heckman2[I, ] <- c(heck2$estimate[3, 1:2], heck2$estimate[4, 1:2])
  heckman3[I, ] <- c(heck3$estimate[4, 1:2], heck3$estimate[5, 1:2])
  heckman4[I, ] <- c(heck4$estimate[4, 1:2], heck4$estimate[5, 1:2])
  
  ## Store diagnostics.
  diagnostics[I, 1] <- glm(R2 ~ 1, family = binomial)$deviance - glm(R2 ~ Z2, family = binomial)$deviance
  diagnostics[I, 2] <- glm(R2 ~ X2, family = binomial)$deviance - glm(R2 ~ X2 + Z2, family = binomial)$deviance
  diagnostics[I, 3] <- glm(R2 ~ X2 + Y2, family = binomial)$deviance - glm(R2 ~ X2 + Y2 + Z2, family = binomial)$deviance
  diagnostics[I, 4] <- mean(R2)
  #diagnostics[I, 5] <- mean(R2[X2 > 0]) - mean(R2[X2 < 0])
  #diagnostics[I, 6] <- (mean(R2[X2 > 0]) * (1 - mean(R2[X2 < 0]))) / (mean(R2[X2 < 0]) * (1 - mean(R2[X2 > 0])))
  diagnostics[I, 7] <- mean(R2[Y2 > 1]) - mean(R2[Y2 < 1])
  diagnostics[I, 8] <- (mean(R2[Y2 > 1]) * (1 - mean(R2[Y2 < 1]))) / (mean(R2[Y2 < 1]) * (1 - mean(R2[Y2 > 1])))
  diagnostics[I, 9] <- ipw.est1$w.min
  diagnostics[I, 10] <- ipw.est1$w.max
  diagnostics[I, 11] <- ipw.est2$w.min
  diagnostics[I, 12] <- ipw.est2$w.max
  diagnostics[I, 13] <- mean(R2[Z2 > 0]) - mean(R2[Z2 < 0])
  diagnostics[I, 14] <- (mean(R2[Z2 > 0]) * (1 - mean(R2[Z2 < 0]))) / (mean(R2[Z2 < 0]) * (1 - mean(R2[Z2 > 0])))
  diagnostics[I, 15] <- summary(lm(X1 ~ G1))$f[1]
  diagnostics[I, 16] <- beta.x^2 * var(G1) / var(X1)
  diagnostics[I, 17] <- beta.x^2 * beta.y^2 * var(G2) / var(Y2)
  diagnostics[I, 18] <- summary(lm(X1 ~ G1))$r.sq
  diagnostics[I, 19] <- summary(lm(Y2 ~ G2))$r.sq
  diagnostics[I, 20] <- 1 - glm(R2 ~ Z2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 21] <- 1 - glm(R2 ~ X2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 22] <- 1 - glm(R2 ~ Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 23] <- 1 - glm(R2 ~ X2 + Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  
  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 100 == 0) save(beta.x, gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n1, n2, gx, oracle, naive,
                          partial, full, diagnostics, ipw1, ipw2, ipw3, 
                          heckman1, heckman2, heckman3, heckman4, 
                          file = filename)
  
}
