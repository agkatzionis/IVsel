
##########   MR SIMULATION 1   ##########

## One Instrument - One Sample.

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load R packages.
library(sampleSelection)
library(survey)

## Results will be saved here.
filename <- "MR1_results.RData"

## Set up the data generation.
seed <- 67199383
iter <- 10000
n <- 10000

## Set values for the main simulation parameters.
beta.x <- sqrt(2/19)   ## 5% genetic instrument strength.
gamma.x <- 1   ## Semi-arbitrary (50% of residual variation in X is due to confounding).
beta.y <- 0.2   ## Tuned to facilitate power comparison.
gamma.y <- 1   ## Semi-arbitrary (50% of residual variation in Y is due to confounding).
alpha.r <- 0   ## Tuned so that 50% of participants are selected.
beta.r <- 0   ## No X-R effect.
gamma.r <- 0.5   ## McFadden's R2 = 0.02.
delta.r <- 1   ## Semi-arbitrary (McFadden's R2 = 0.20-0.25).

## Make sure the proportions of genetic variation 
## in the risk factor and outcome are reasonable.
gv.x <- (beta.x^2 * 1) / (beta.x^2 + gamma.x^2 + 1)
gv.y <- (beta.y^2 * beta.x^2) / (beta.y^2 * beta.x^2 + (beta.y * gamma.x + gamma.y)^2 + beta.y^2 + 1)

## Store G-X associations here.
gx <- matrix(0, iter, 4); colnames(gx) <- c("Intercept", "se.int", "Beta", "se.beta")

## Store G-Y associations here.
oracle <- matrix(0, iter, 4); colnames(oracle) <- c("Intercept", "se.int", "Beta", "se.beta")
naive <- oracle; ipw1 <- oracle; ipw2 <- oracle; ipw3 <- oracle; ipw4 <- oracle; ipw5 <- oracle
heckman1 <- oracle; heckman2 <- oracle; heckman3 <- oracle
heckman4 <- oracle; heckman5 <- oracle; heckman6 <- oracle
partial <- oracle; full <- oracle
svyipw1 <- oracle; svyipw2 <- oracle; svyipw3 <- oracle; svyipw4 <- oracle; svyipw5 <- oracle

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
  gx.est <- summary(lm(X ~ G))
  
  ## Naive and oracle estimates.
  naive.est <- summary(lm(Y[R == 1] ~ G[R == 1]))
  oracle.est <- summary(lm(Y ~ G))
  
  ## Inverse probability weighting.
  ipw.est1 <- ipw.linear(X = G, Y = Y, R = R)
  ipw.est2 <- ipw.linear(X = G, Y = Y, R = R, Z = Z)
  ipw.est3 <- ipw.linear(X = G, Y = Y, R = R, Z = Z, drop.x = TRUE)
  ipw.est4 <- ipw.linear(X = G, Y = Y, R = R, Z = X)
  ipw.est5 <- ipw.linear(X = G, Y = Y, R = R, Z = cbind(X, Z))
  
  ## Heckman's method.
  Yobs <- Y;  Yobs[R == 0] <- NA
  heck1 <- summary(selection(R ~ Z, Yobs ~ G, method = "2step"))
  heck2 <- summary(selection(R ~ Z, Yobs ~ G, method = "ml"))
  heck3 <- summary(selection(R ~ G + Z, Yobs ~ G, method = "2step"))
  heck4 <- summary(selection(R ~ G + Z, Yobs ~ G, method = "ml"))
  heck5 <- summary(selection(R ~ G + X + Z, Yobs ~ G, method = "2step"))
  heck6 <- summary(selection(R ~ G + X + Z, Yobs ~ G, method = "ml"))
  
  ## Partial likelihood MLE.
  partial.opt1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z)
  partial.opt2 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = Y, alpha.hat = partial.opt1$par, hessian = TRUE)
  partial.sd <- sqrt(diag(solve(- partial.opt2$hessian))[1:2])
  
  ## Full likelihood MLE.
  full.opt <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = Y, hessian = TRUE)
  full.sd <- sqrt(diag(solve(- full.opt$hessian))[1:2])
  
  ## Store G-X results.
  gx[I, ] <- c(gx.est$coef[1, 1:2], gx.est$coef[2, 1:2])
  
  ## Store G-Y results.
  oracle[I, ] <- c(oracle.est$coef[1, 1], oracle.est$coef[1, 2], oracle.est$coef[2, 1], oracle.est$coef[2, 2])
  naive[I, ] <- c(naive.est$coef[1, 1], naive.est$coef[1, 2], naive.est$coef[2, 1], naive.est$coef[2, 2])
  ipw1[I, ] <- c(ipw.est1$est[1, 1], ipw.est1$est[1, 2], ipw.est1$est[2, 1], ipw.est1$est[2, 2])
  ipw2[I, ] <- c(ipw.est2$est[1, 1], ipw.est2$est[1, 2], ipw.est2$est[2, 1], ipw.est2$est[2, 2])
  ipw3[I, ] <- c(ipw.est3$est[1, 1], ipw.est3$est[1, 2], ipw.est3$est[2, 1], ipw.est3$est[2, 2])
  ipw4[I, ] <- c(ipw.est4$est[1, 1], ipw.est4$est[1, 2], ipw.est4$est[2, 1], ipw.est4$est[2, 2])
  ipw5[I, ] <- c(ipw.est5$est[1, 1], ipw.est5$est[1, 2], ipw.est5$est[2, 1], ipw.est5$est[2, 2])
  partial[I, ] <- c(partial.opt2$par[1], partial.sd[1], partial.opt2$par[2], partial.sd[2])
  full[I, ] <- c(full.opt$par[1], full.sd[1], full.opt$par[2], full.sd[2])

  ## Store more results.
  heckman1[I, ] <- c(heck1$estimate[3, 1:2], heck1$estimate[4, 1:2])
  heckman2[I, ] <- c(heck2$estimate[3, 1:2], heck2$estimate[4, 1:2])
  heckman3[I, ] <- c(heck3$estimate[4, 1:2], heck3$estimate[5, 1:2])
  heckman4[I, ] <- c(heck4$estimate[4, 1:2], heck4$estimate[5, 1:2])
  heckman5[I, ] <- c(heck5$estimate[5, 1:2], heck5$estimate[6, 1:2])
  heckman6[I, ] <- c(heck6$estimate[5, 1:2], heck6$estimate[6, 1:2])
  
  ## Store IPW results with robust standard errors.
  svyipw1[I, ] <- c(ipw.est1$svyest[1, 1], ipw.est1$svyest[1, 2], ipw.est1$svyest[2, 1], ipw.est1$svyest[2, 2])
  svyipw2[I, ] <- c(ipw.est2$svyest[1, 1], ipw.est2$svyest[1, 2], ipw.est2$svyest[2, 1], ipw.est2$svyest[2, 2])
  svyipw3[I, ] <- c(ipw.est3$svyest[1, 1], ipw.est3$svyest[1, 2], ipw.est3$svyest[2, 1], ipw.est3$svyest[2, 2])
  svyipw4[I, ] <- c(ipw.est4$svyest[1, 1], ipw.est4$svyest[1, 2], ipw.est4$svyest[2, 1], ipw.est4$svyest[2, 2])
  svyipw5[I, ] <- c(ipw.est5$svyest[1, 1], ipw.est5$svyest[1, 2], ipw.est5$svyest[2, 1], ipw.est5$svyest[2, 2])
  
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
  if (I %% 100 == 0) save(beta.x, gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n, gx, 
                          oracle, naive, partial, full, ipw1, ipw2, ipw3, ipw4, ipw5, 
                          heckman1, heckman2, heckman3, heckman4, heckman5, heckman6, diagnostics,
                          svyipw1, svyipw2, svyipw3, svyipw4, svyipw5, file = filename)
  
}
