
##########   MR SIMULATION 12   ##########

## One Instrument - Two Samples.

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)

## Results will be saved here.
filename <- "MR12_results.RData"

## Set up the data generation.
seed <- 67210383
iter <- 10000
n1 <- 10000
n2 <- 10000

## Set values for the main simulation parameters.
beta.x <- sqrt(2 / 19)
gamma.x <- 1
beta.y <- 0   ## Null causal effect.
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
diagnostics <- matrix(0, iter, 29)
colnames(diagnostics) <- c("Strength 1", "Strength 2", "Strength | X", "Strength | X, Y", "Prop Sel X", "Prop Sel Y", "MD X", "OR X", "MD Y", "OR Y", 
                           "Wmin | X", "Wmax | X", "Wmin | X, Z", "Wmax | X, Z", "MD Z1", "OR Z1", "MD Z2", "OR Z2", "F-stat G", "GV% X", "GV% Y",
                           "R2 (X ~ G)", "R2 (Y ~ G)", "McFadden Z1", "McFadden Z2", "McFadden X", "McFadden Y", "McFadden XY1", "McFadden XY2")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate the first sample.
  U1 <- rnorm(n1, 0, 1)
  Z1 <- rnorm(n1, 0, 1)
  G1 <- rnorm(n1, 0, 1)
  X1 <- beta.x * G1 + gamma.x * U1 + rnorm(n1, 0, 1)
  Y1 <- beta.y * X1 + gamma.y * U1 + rnorm(n1, 0, 1)
  
  ## Simulate the second sample.
  U2 <- rnorm(n2, 0, 1)
  Z2 <- rnorm(n2, 0, 1)
  G2 <- rnorm(n2, 0, 1)
  X2 <- beta.x * G2 + gamma.x * U2 + rnorm(n2, 0, 1)
  Y2 <- beta.y * X2 + gamma.y * U2 + rnorm(n2, 0, 1)
  
  ## Simulate the selection coefficient.
  R1.probs <- expit(alpha.r + beta.r * X1 + gamma.r * Z1 + delta.r * Y1)
  R1 <- rbinom(n1, 1, R1.probs)
  R2.probs <- expit(alpha.r + beta.r * X2 + gamma.r * Z2 + delta.r * Y2)
  R2 <- rbinom(n2, 1, R2.probs)
  
  ## G-X association estimates.
  
  ## Naive and oracle estimates.
  naive.est.x <- summary(lm(X1[R1 == 1] ~ G1[R1 == 1]))
  oracle.est.x <- summary(lm(X1 ~ G1))
  
  ## Inverse probability weighting.
  ipw.est1.x <- ipw.linear(X = G1, Y = X1, R = R1)
  ipw.est2.x <- ipw.linear(X = G1, Y = X1, R = R1, Z = Z1)
  ipw.est3.x <- ipw.linear(X = G1, Y = X1, R = R1, Z = Z1, drop.x = TRUE)
  
  ## Heckman's method.
  X1obs <- X1;  X1obs[R1 == 0] <- NA
  heck1.x <- summary(selection(R1 ~ Z1, X1obs ~ G1, method = "2step"))
  heck2.x <- summary(selection(R1 ~ Z1, X1obs ~ G1, method = "ml"))
  heck3.x <- summary(selection(R1 ~ G1 + Z1, X1obs ~ G1, method = "2step"))
  heck4.x <- summary(selection(R1 ~ G1 + Z1, X1obs ~ G1, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1.x <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G1, R = R1, Z = Z1)
  partial.opt2.x <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G1, R = R1, Z = Z1, Y = X1, alpha.hat = partial.opt1.x$par, hessian = TRUE)
  partial.sd.x <- sqrt(diag(solve(- partial.opt2.x$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt.x <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G1, R = R1, Z = Z1, Y = X1, hessian = TRUE)
  full.sd.x <- sqrt(diag(solve(- full.opt.x$hessian))[1:2])
  
  ## G-X association estimates.
  
  ## Naive and oracle estimates.
  naive.est.y <- summary(lm(Y2[R2 == 1] ~ G2[R2 == 1]))
  oracle.est.y <- summary(lm(Y2 ~ G2))
  
  ## Inverse probability weighting.
  ipw.est1.y <- ipw.linear(X = G2, Y = Y2, R = R2)
  ipw.est2.y <- ipw.linear(X = G2, Y = Y2, R = R2, Z = Z2)
  ipw.est3.y <- ipw.linear(X = G2, Y = Y2, R = R2, Z = Z2, drop.x = TRUE)
  
  ## Heckman's method.
  Y2obs <- Y2;  Y2obs[R2 == 0] <- NA
  heck1.y <- summary(selection(R2 ~ Z2, Y2obs ~ G2, method = "2step"))
  heck2.y <- summary(selection(R2 ~ Z2, Y2obs ~ G2, method = "ml"))
  heck3.y <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2, method = "2step"))
  heck4.y <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1.y <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2)
  partial.opt2.y <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2, Y = Y2, alpha.hat = partial.opt1.y$par, hessian = TRUE)
  partial.sd.y <- sqrt(diag(solve(- partial.opt2.y$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt.y <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G2, R = R2, Z = Z2, Y = Y2, hessian = TRUE)
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
  diagnostics[I, 1] <- glm(R1 ~ 1, family = binomial)$deviance - glm(R1 ~ Z1, family = binomial)$deviance
  diagnostics[I, 2] <- glm(R2 ~ 1, family = binomial)$deviance - glm(R2 ~ Z2, family = binomial)$deviance
  diagnostics[I, 3] <- glm(R1 ~ X1, family = binomial)$deviance - glm(R1 ~ X1 + Z1, family = binomial)$deviance
  diagnostics[I, 4] <- glm(R2 ~ X2 + Y2, family = binomial)$deviance - glm(R2 ~ X2 + Y2 + Z2, family = binomial)$deviance
  diagnostics[I, 5] <- mean(R1)
  diagnostics[I, 6] <- mean(R2)
  diagnostics[I, 7] <- mean(R1[X1 > 0]) - mean(R1[X1 < 0])
  diagnostics[I, 8] <- (mean(R1[X1 > 0]) * (1 - mean(R1[X1 < 0]))) / (mean(R1[X1 < 0]) * (1 - mean(R1[X1 > 0])))
  diagnostics[I, 9] <- mean(R2[Y2 > 1]) - mean(R2[Y2 < 1])
  diagnostics[I, 10] <- (mean(R2[Y2 > 1]) * (1 - mean(R2[Y2 < 1]))) / (mean(R2[Y2 < 1]) * (1 - mean(R2[Y2 > 1])))
  diagnostics[I, 11] <- ipw.est1.y$w.min
  diagnostics[I, 12] <- ipw.est1.y$w.max
  diagnostics[I, 13] <- ipw.est2.x$w.min
  diagnostics[I, 14] <- ipw.est2.x$w.max
  diagnostics[I, 15] <- mean(R1[Z1 > 0]) - mean(R1[Z1 < 0])
  diagnostics[I, 16] <- (mean(R1[Z1 > 0]) * (1 - mean(R1[Z1 < 0]))) / (mean(R1[Z1 < 0]) * (1 - mean(R1[Z1 > 0])))
  diagnostics[I, 17] <- mean(R2[Z2 > 0]) - mean(R2[Z2 < 0])
  diagnostics[I, 18] <- (mean(R2[Z2 > 0]) * (1 - mean(R2[Z2 < 0]))) / (mean(R2[Z2 < 0]) * (1 - mean(R2[Z2 > 0])))
  diagnostics[I, 19] <- summary(lm(X1 ~ G1))$f[1]
  diagnostics[I, 20] <- beta.x^2 * var(G1) / var(X1)
  diagnostics[I, 21] <- beta.x^2 * beta.y^2 * var(G2) / var(Y2)
  diagnostics[I, 22] <- summary(lm(X1 ~ G1))$r.sq
  diagnostics[I, 23] <- summary(lm(Y2 ~ G2))$r.sq
  diagnostics[I, 24] <- 1 - glm(R1 ~ Z1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 25] <- 1 - glm(R2 ~ Z2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 26] <- 1 - glm(R1 ~ X1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 27] <- 1 - glm(R2 ~ Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 28] <- 1 - glm(R1 ~ X1 + Y1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 29] <- 1 - glm(R2 ~ X2 + Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  
  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 100 == 0) save(beta.x, gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n1, n2, diagnostics, 
                          oracle.x, naive.x, partial.x, full.x, ipw1.x, ipw2.x, ipw3.x, 
                          heckman1.x, heckman2.x, heckman3.x, heckman4.x, 
                          oracle.y, naive.y, partial.y, full.y, ipw1.y, ipw2.y, ipw3.y, 
                          heckman1.y, heckman2.y, heckman3.y, heckman4.y, 
                          file = filename)
  
}
