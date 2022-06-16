
##########   MR SIMULATION 21   ##########

## Many Instruments - One Sample.

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)
library(ivreg)
library(MendelianRandomization)
library(truncnorm)

## Results will be saved here.
filename <- "MR21_results.RData"

## Set up the data generation.
seed <- 67219383
iter <- 1000
n1 <- 10000
n2 <- 10000
K <- 10

## Set values for the main simulation parameters.
gamma.x <- 1
beta.y <- 0.2
gamma.y <- 1
beta.r <- 0.5   ## Semi-arbitrary (joint McFadden's R2 = 0.20-0.25).
gamma.r <- 0.5
delta.r <- 0.5   ## Semi-arbitrary (joint McFadden's R2 = 0.20-0.25).
alpha.r <- - (beta.r + beta.y * delta.r) / 2 * dnorm(3) / (1 - pnorm(3))

## Beta.x and gv will be generated within the iteration loop.
## Proportions of genetic variation depend on allele frequencies.

## Store betas and effect allele frequencies here.
eafs <- matrix(0, iter, K)
betas <- matrix(0, iter, K)

## Store G-X summary statistics here.
bx.oracle.all <- matrix(0, iter, K); sx.oracle.all <- matrix(0, iter, K)
bx.naive.all <- matrix(0, iter, K); sx.naive.all <- matrix(0, iter, K)
bx.ipw1.all <- matrix(0, iter, K); sx.ipw1.all <- matrix(0, iter, K)
bx.ipw2.all <- matrix(0, iter, K); sx.ipw2.all <- matrix(0, iter, K)
bx.ipw3.all <- matrix(0, iter, K); sx.ipw3.all <- matrix(0, iter, K)
bx.heckman1.all <- matrix(0, iter, K); sx.heckman1.all <- matrix(0, iter, K)
bx.heckman2.all <- matrix(0, iter, K); sx.heckman2.all <- matrix(0, iter, K)
bx.heckman3.all <- matrix(0, iter, K); sx.heckman3.all <- matrix(0, iter, K)
bx.heckman4.all <- matrix(0, iter, K); sx.heckman4.all <- matrix(0, iter, K)
bx.heckman5.all <- matrix(0, iter, K); sx.heckman5.all <- matrix(0, iter, K)
bx.heckman6.all <- matrix(0, iter, K); sx.heckman6.all <- matrix(0, iter, K)
bx.partial1.all <- matrix(0, iter, K); sx.partial1.all <- matrix(0, iter, K)
bx.partial2.all <- matrix(0, iter, K); sx.partial2.all <- matrix(0, iter, K)
bx.partial3.all <- matrix(0, iter, K); sx.partial3.all <- matrix(0, iter, K)
bx.full.all <- matrix(0, iter, K); sx.full.all <- matrix(0, iter, K)

## Store G-Y summary statistics here.
by.oracle.all <- matrix(0, iter, K); sy.oracle.all <- matrix(0, iter, K)
by.naive.all <- matrix(0, iter, K); sy.naive.all <- matrix(0, iter, K)
by.ipw1.all <- matrix(0, iter, K); sy.ipw1.all <- matrix(0, iter, K)
by.ipw2.all <- matrix(0, iter, K); sy.ipw2.all <- matrix(0, iter, K)
by.ipw3.all <- matrix(0, iter, K); sy.ipw3.all <- matrix(0, iter, K)
by.heckman1.all <- matrix(0, iter, K); sy.heckman1.all <- matrix(0, iter, K)
by.heckman2.all <- matrix(0, iter, K); sy.heckman2.all <- matrix(0, iter, K)
by.heckman3.all <- matrix(0, iter, K); sy.heckman3.all <- matrix(0, iter, K)
by.heckman4.all <- matrix(0, iter, K); sy.heckman4.all <- matrix(0, iter, K)
by.heckman5.all <- matrix(0, iter, K); sy.heckman5.all <- matrix(0, iter, K)
by.heckman6.all <- matrix(0, iter, K); sy.heckman6.all <- matrix(0, iter, K)
by.partial1.all <- matrix(0, iter, K); sy.partial1.all <- matrix(0, iter, K)
by.partial2.all <- matrix(0, iter, K); sy.partial2.all <- matrix(0, iter, K)
by.partial3.all <- matrix(0, iter, K); sy.partial3.all <- matrix(0, iter, K)
by.full.all <- matrix(0, iter, K); sy.full.all <- matrix(0, iter, K)

## Store IVW MR estimates here.
oracle.s <- matrix(0, iter, 2); colnames(oracle.s) <- c("Causal", "se.causal")
naive.s <- oracle.s; ipw1.s <- oracle.s; ipw2.s <- oracle.s
ipw3.s <- oracle.s ; full.s <- oracle.s; partial1.s <- oracle.s
partial2.s <- oracle.s; partial3.s <- oracle.s; heckman1.s <- oracle.s
heckman2.s <- oracle.s; heckman3.s <- oracle.s; heckman4.s <- oracle.s
heckman5.s <- oracle.s; heckman6.s <- oracle.s

## Store additional diagnostics here.
diagnostics <- matrix(0, iter, 25)
colnames(diagnostics) <- c("Strength 1", "Strength 2", "Strength | X", "Strength | X, Y", "Prop Sel X", "Prop Sel Y", "MD X", "OR X", "MD Y", "OR Y", 
                           "MD Z1", "OR Z1", "MD Z2", "OR Z2", "F-stat G", "GV% X", "GV% Y", "R2 (X ~ G)", "R2 (Y ~ G)", 
                           "McFadden Z1", "McFadden Z2", "McFadden X", "McFadden Y", "McFadden XY1", "McFadden XY2")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate betas and effect allele frequencies.
  beta.x <- rtruncnorm(K, a = 0.15, mean = 0, sd = 0.05)
  eaf <- runif(K, 0.1, 0.9)
  
  ## Generate the first sample.
  U1 <- rnorm(n1, 0, 1)
  Z1 <- rnorm(n1, 0, 1)
  G1 <- matrix(0, n1, K)
  for (j in 1:K) G1[, j] <- rbinom(n1, 2, eaf[j])
  X1 <- as.vector(G1 %*% beta.x) + gamma.x * U1 + rnorm(n1, 0, 1)
  Y1 <- beta.y * X1 + gamma.y * U1 + rnorm(n1, 0, 1)
  
  ## Simulate the first sample's selection coefficient.
  R1.probs <- expit(alpha.r + beta.r * X1 + gamma.r * Z1 + delta.r * Y1)
  R1 <- rbinom(n1, 1, R1.probs)
  X1obs <- X1; X1obs[R1 == 0] <- NA
  
  ## Generate the second sample.
  U2 <- rnorm(n2, 0, 1)
  Z2 <- rnorm(n2, 0, 1)
  G2 <- matrix(0, n2, K)
  for (j in 1:K) G2[, j] <- rbinom(n2, 2, eaf[j])
  X2 <- as.vector(G2 %*% beta.x) + gamma.x * U2 + rnorm(n2, 0, 1)
  Y2 <- beta.y * X2 + gamma.y * U2 + rnorm(n2, 0, 1)
  
  ## Simulate the second sample's selection coefficient.
  R2.probs <- expit(alpha.r + beta.r * X2 + gamma.r * Z2 + delta.r * Y2)
  R2 <- rbinom(n2, 1, R2.probs)
  Y2obs <- Y2; Y2obs[R2 == 0] <- NA
  
  ## Summary statistics and IVW estimates.
  
  ## Oracle summary statistics.
  bx.oracle <- rep(0, K); sx.oracle <- rep(0, K)
  by.oracle <- rep(0, K); sy.oracle <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(lm(X1 ~ G1[, i]))$coefficients
    bx.oracle[i] <- fit1[2, 1]
    sx.oracle[i] <- fit1[2, 2]
    fit2 <- summary(lm(Y2 ~ G2[, i]))$coefficients
    by.oracle[i] <- fit2[2, 1]
    sy.oracle[i] <- fit2[2, 2]
  }
  
  ## Oracle causal effect estimate.
  oracle.est <- mr_ivw(mr_input(bx = bx.oracle, bxse = sx.oracle, by = by.oracle, byse = sy.oracle))
  
  ## Complete-case summary statistics.
  bx.naive <- rep(0, K); sx.naive <- rep(0, K)
  by.naive <- rep(0, K); sy.naive <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(lm(X1[R1 == 1] ~ G1[R1 == 1, i]))$coefficients
    bx.naive[i] <- fit1[2, 1]
    sx.naive[i] <- fit1[2, 2]
    fit2 <- summary(lm(Y2[R2 == 1] ~ G2[R2 == 1, i]))$coefficients
    by.naive[i] <- fit2[2, 1]
    sy.naive[i] <- fit2[2, 2]
  }
  
  ## Complete-case causal effect estimate.
  naive.est <- mr_ivw(mr_input(bx = bx.naive, bxse = sx.naive, by = by.naive, byse = sy.naive))
  
  ## IPW (G) summary statistics.
  bx.ipw1 <- rep(0, K); sx.ipw1 <- rep(0, K)
  by.ipw1 <- rep(0, K); sy.ipw1 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G1[, i], Y = X1, R = R1, Z = G1[, which(!(1:K %in% i))])
    bx.ipw1[i] <- fit1$est[2, 1]
    sx.ipw1[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G2[, i], Y = Y2, R = R2, Z = G2[, which(!(1:K %in% i))])
    by.ipw1[i] <- fit2$est[2, 1]
    sy.ipw1[i] <- fit2$est[2, 2]
  }
  
  ## IPW (G) causal effect estimate.
  ipw.est1 <- mr_ivw(mr_input(bx = bx.ipw1, bxse = sx.ipw1, by = by.ipw1, byse = sy.ipw1))
  
  ## IPW (G, Z) summary statistics.
  bx.ipw2 <- rep(0, K); sx.ipw2 <- rep(0, K)
  by.ipw2 <- rep(0, K); sy.ipw2 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G1[, i], Y = X1, R = R1, Z = cbind(Z1, G1[, which(!(1:K %in% i))]))
    bx.ipw2[i] <- fit1$est[2, 1]
    sx.ipw2[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G2[, i], Y = Y2, R = R2, Z = cbind(Z2, G2[, which(!(1:K %in% i))]))
    by.ipw2[i] <- fit2$est[2, 1]
    sy.ipw2[i] <- fit2$est[2, 2]
  }
  
  ## IPW (G, Z) causal effect estimate.
  ipw.est2 <- mr_ivw(mr_input(bx = bx.ipw2, bxse = sx.ipw2, by = by.ipw2, byse = sy.ipw2))
  
  ## IPW (Z) summary statistics.
  bx.ipw3 <- rep(0, K); sx.ipw3 <- rep(0, K)
  by.ipw3 <- rep(0, K); sy.ipw3 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G1[, i], Y = X1, R = R1, Z = Z1, drop.x = TRUE)
    bx.ipw3[i] <- fit1$est[2, 1]
    sx.ipw3[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G2[, i], Y = Y2, R = R2, Z = Z2, drop.x = TRUE)
    by.ipw3[i] <- fit2$est[2, 1]
    sy.ipw3[i] <- fit2$est[2, 2]
  }
  
  ## IPW (Z) causal effect estimate.
  ipw.est3 <- mr_ivw(mr_input(bx = bx.ipw3, bxse = sx.ipw3, by = by.ipw3, byse = sy.ipw3))
  
  ## Heckman-2SLS (Z) summary statistics.
  bx.heck1 <- rep(0, K); sx.heck1 <- rep(0, K)
  by.heck1 <- rep(0, K); sy.heck1 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ Z1, X1obs ~ G1[, i], method = "2step"))[[2]]
    bx.heck1[i] <- fit1[4, 1]
    sx.heck1[i] <- fit1[4, 2]
    fit2 <- summary(selection(R2 ~ Z2, Y2obs ~ G2[, i], method = "2step"))[[2]]
    by.heck1[i] <- fit2[4, 1]
    sy.heck1[i] <- fit2[4, 2]
  }
  
  ## Heckman-2SLS (Z) causal effect estimate.
  heck.est1 <- mr_ivw(mr_input(bx = bx.heck1, bxse = sx.heck1, by = by.heck1, byse = sy.heck1))
  
  ## Heckman-ML (Z) summary statistics.
  bx.heck2 <- rep(0, K); sx.heck2 <- rep(0, K)
  by.heck2 <- rep(0, K); sy.heck2 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ Z1, X1obs ~ G1[, i], method = "ml"))$estimate
    bx.heck2[i] <- fit1[4, 1]
    sx.heck2[i] <- fit1[4, 2]
    fit2 <- summary(selection(R2 ~ Z2, Y2obs ~ G2[, i], method = "ml"))$estimate
    by.heck2[i] <- fit2[4, 1]
    sy.heck2[i] <- fit2[4, 2]
  }
  
  ## Heckman-ML (Z) causal effect estimate.
  heck.est2 <- mr_ivw(mr_input(bx = bx.heck2, bxse = sx.heck2, by = by.heck2, byse = sy.heck2))
  
  ## Heckman-2SLS (G, Z) summary statistics.
  bx.heck3 <- rep(0, K); sx.heck3 <- rep(0, K)
  by.heck3 <- rep(0, K); sy.heck3 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ G1 + Z1, X1obs ~ G1[, i], method = "2step"))[[2]]
    bx.heck3[i] <- fit1[K + 4, 1]
    sx.heck3[i] <- fit1[K + 4, 2]
    fit2 <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2[, i], method = "2step"))[[2]]
    by.heck3[i] <- fit2[K + 4, 1]
    sy.heck3[i] <- fit2[K + 4, 2]
  }
  
  ## Heckman-2SLS (G, Z) causal effect estimate.
  heck.est3 <- mr_ivw(mr_input(bx = bx.heck3, bxse = sx.heck3, by = by.heck3, byse = sy.heck3))
  
  ## Heckman-ML (G, Z) summary statistics.
  bx.heck4 <- rep(0, K); sx.heck4 <- rep(0, K)
  by.heck4 <- rep(0, K); sy.heck4 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ G1 + Z1, X1obs ~ G1[, i], method = "ml"))$estimate
    bx.heck4[i] <- fit1[K + 4, 1]
    sx.heck4[i] <- fit1[K + 4, 2]
    fit2 <- summary(selection(R2 ~ G2 + Z2, Y2obs ~ G2[, i], method = "ml"))$estimate
    by.heck4[i] <- fit2[K + 4, 1]
    sy.heck4[i] <- fit2[K + 4, 2]
  }
  
  ## Heckman-ML (G, Z) causal effect estimate.
  heck.est4 <- mr_ivw(mr_input(bx = bx.heck4, bxse = sx.heck4, by = by.heck4, byse = sy.heck4))
  
  ## Heckman-2SLS (G, Z) summary statistics.
  bx.heck5 <- rep(0, K); sx.heck5 <- rep(0, K)
  by.heck5 <- rep(0, K); sy.heck5 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ G1[, i] + Z1, X1obs ~ G1[, i], method = "2step"))[[2]]
    bx.heck5[i] <- fit1[5, 1]
    sx.heck5[i] <- fit1[5, 2]
    fit2 <- summary(selection(R2 ~ G2[, i] + Z2, Y2obs ~ G2[, i], method = "2step"))[[2]]
    by.heck5[i] <- fit2[5, 1]
    sy.heck5[i] <- fit2[5, 2]
  }
  
  ## Heckman-2SLS (G, Z) causal effect estimate.
  heck.est5 <- mr_ivw(mr_input(bx = bx.heck5, bxse = sx.heck5, by = by.heck5, byse = sy.heck5))
  
  ## Heckman-ML (G, Z) summary statistics.
  bx.heck6 <- rep(0, K); sx.heck6 <- rep(0, K)
  by.heck6 <- rep(0, K); sy.heck6 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R1 ~ G1[, i] + Z1, X1obs ~ G1[, i], method = "ml"))$estimate
    bx.heck6[i] <- fit1[5, 1]
    sx.heck6[i] <- fit1[5, 2]
    fit2 <- summary(selection(R2 ~ G2[, i] + Z2, Y2obs ~ G2[, i], method = "ml"))$estimate
    by.heck6[i] <- fit2[5, 1]
    sy.heck6[i] <- fit2[5, 2]
  }
  
  ## Heckman-ML (G, Z) causal effect estimate.
  heck.est6 <- mr_ivw(mr_input(bx = bx.heck6, bxse = sx.heck6, by = by.heck6, byse = sy.heck6))
  
  ## TTW-partial (default) summary statistics.
  bx.partial1 <- rep(0, K); sx.partial1 <- rep(0, K)
  by.partial1 <- rep(0, K); sy.partial1 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1.1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = Z1)
    partial.opt2 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = Z1, Y = X1, alpha.hat = partial.opt1.1$par, hessian = TRUE)
    bx.partial1[i] <- partial.opt2$par[2]
    sx.partial1[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt1.2 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = Z2)
    partial.opt3 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = Z2, Y = Y2, alpha.hat = partial.opt1.2$par, hessian = TRUE)
    by.partial1[i] <- partial.opt3$par[2]
    sy.partial1[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (default) causal effect estimate.
  partial.est1 <- mr_ivw(mr_input(bx = bx.partial1, bxse = sx.partial1, by = by.partial1, byse = sy.partial1))
  
  ## TTW-full (default) summary statistics.
  bx.full <- rep(0, K); sx.full <- rep(0, K)
  by.full <- rep(0, K); sy.full <- rep(0, K)
  for (i in 1:K) {
    full.opt1 <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = Z1, Y = X1, hessian = TRUE)
    bx.full[i] <- full.opt1$par[2]
    sx.full[i] <- sqrt(diag(solve(- full.opt1$hessian))[2])
    full.opt2 <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = Z2, Y = Y2, hessian = TRUE)
    by.full[i] <- full.opt2$par[2]
    sy.full[i] <- sqrt(diag(solve(- full.opt2$hessian))[2])
  }
  
  ## TTW-full (default) causal effect estimate.
  full.est <- mr_ivw(mr_input(bx = bx.full, bxse = sx.full, by = by.full, byse = sy.full))
  
  ## TTW-partial (G-prop) summary statistics.
  bx.partial2 <- rep(0, K); sx.partial2 <- rep(0, K)
  by.partial2 <- rep(0, K); sy.partial2 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1.1 <- optim(par = rep(0, 12), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = cbind(Z1, G1[, which(!(1:K %in% i))]))
    partial.opt2 <- optim(par = rep(0, 5), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = cbind(Z1, G1[, which(!(1:K %in% i))]), Y = X1, alpha.hat = partial.opt1.1$par, hessian = TRUE)
    bx.partial2[i] <- partial.opt2$par[2]
    sx.partial2[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt1.2 <- optim(par = rep(0, 12), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = cbind(Z2, G2[, which(!(1:K %in% i))]))
    partial.opt3 <- optim(par = rep(0, 5), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = cbind(Z2, G2[, which(!(1:K %in% i))]), Y = Y2, alpha.hat = partial.opt1.2$par, hessian = TRUE)
    by.partial2[i] <- partial.opt3$par[2]
    sy.partial2[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (G-prop) causal effect estimate.
  partial.est2 <- mr_ivw(mr_input(bx = bx.partial2, bxse = sx.partial2, by = by.partial2, byse = sy.partial2))
  
  ## TTW-partial (G-sel) summary statistics.
  bx.partial3 <- rep(0, K); sx.partial3 <- rep(0, K)
  by.partial3 <- rep(0, K); sy.partial3 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1.1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = Z1)
    partial.opt2 <- optim(par = rep(0, 14), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G1[, i], R = R1, Z = Z1, Y = X1, add.sel = G1[, which(!(1:K %in% i))], alpha.hat = partial.opt1.1$par, hessian = TRUE)
    bx.partial3[i] <- partial.opt2$par[2]
    sx.partial3[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt1.2 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = Z2)
    partial.opt3 <- optim(par = rep(0, 14), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G2[, i], R = R2, Z = Z2, Y = Y2, add.sel = G2[, which(!(1:K %in% i))], alpha.hat = partial.opt1.2$par, hessian = TRUE)
    by.partial3[i] <- partial.opt3$par[2]
    sy.partial3[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (G-sel) causal effect estimate.
  partial.est3 <- mr_ivw(mr_input(bx = bx.partial3, bxse = sx.partial3, by = by.partial3, byse = sy.partial3))
  
  ## Store betas and allele frequencies.
  eafs[I, ] <- eaf
  betas[I, ] <- beta.x
  
  ## Store G-X summary statistics here.
  bx.oracle.all[I, ] <- bx.oracle; sx.oracle.all[I, ] <- sx.oracle
  bx.naive.all[I, ] <- bx.naive; sx.naive.all[I, ] <- sx.naive
  bx.ipw1.all[I, ] <- bx.ipw1; sx.ipw1.all[I, ] <- sx.ipw1
  bx.ipw2.all[I, ] <- bx.ipw2; sx.ipw2.all[I, ] <- sx.ipw2
  bx.ipw3.all[I, ] <- bx.ipw3; sx.ipw3.all[I, ] <- sx.ipw3
  bx.heckman1.all[I, ] <- bx.heck1; sx.heckman1.all[I, ] <- sx.heck1
  bx.heckman2.all[I, ] <- bx.heck2; sx.heckman2.all[I, ] <- sx.heck2
  bx.heckman3.all[I, ] <- bx.heck3; sx.heckman3.all[I, ] <- sx.heck3
  bx.heckman4.all[I, ] <- bx.heck4; sx.heckman4.all[I, ] <- sx.heck4
  bx.heckman5.all[I, ] <- bx.heck5; sx.heckman5.all[I, ] <- sx.heck5
  bx.heckman6.all[I, ] <- bx.heck6; sx.heckman6.all[I, ] <- sx.heck6
  bx.partial1.all[I, ] <- bx.partial1; sx.partial1.all[I, ] <- sx.partial1
  bx.full.all[I, ] <- bx.full; sx.full.all[I, ] <- sx.full
  bx.partial2.all[I, ] <- bx.partial2; sx.partial2.all[I, ] <- sx.partial2
  bx.partial3.all[I, ] <- bx.partial3; sx.partial3.all[I, ] <- sx.partial3
  
  ## Store G-Y summary statistics.
  by.oracle.all[I, ] <- by.oracle; sy.oracle.all[I, ] <- sy.oracle
  by.naive.all[I, ] <- by.naive; sy.naive.all[I, ] <- sy.naive
  by.ipw1.all[I, ] <- by.ipw1; sy.ipw1.all[I, ] <- sy.ipw1
  by.ipw2.all[I, ] <- by.ipw2; sy.ipw2.all[I, ] <- sy.ipw2
  by.ipw3.all[I, ] <- by.ipw3; sy.ipw3.all[I, ] <- sy.ipw3
  by.heckman1.all[I, ] <- by.heck1; sy.heckman1.all[I, ] <- sy.heck1
  by.heckman2.all[I, ] <- by.heck2; sy.heckman2.all[I, ] <- sy.heck2
  by.heckman3.all[I, ] <- by.heck3; sy.heckman3.all[I, ] <- sy.heck3
  by.heckman4.all[I, ] <- by.heck4; sy.heckman4.all[I, ] <- sy.heck4
  by.heckman5.all[I, ] <- by.heck5; sy.heckman5.all[I, ] <- sy.heck5
  by.heckman6.all[I, ] <- by.heck6; sy.heckman6.all[I, ] <- sy.heck6
  by.partial1.all[I, ] <- by.partial1; sy.partial1.all[I, ] <- sy.partial1
  by.full.all[I, ] <- by.full; sy.full.all[I, ] <- sy.full
  by.partial2.all[I, ] <- by.partial2; sy.partial2.all[I, ] <- sy.partial2
  by.partial3.all[I, ] <- by.partial3; sy.partial3.all[I, ] <- sy.partial3
  
  ## Store causal effect estimates.
  oracle.s[I, ] <- c(oracle.est$Estimate, oracle.est$StdError)
  naive.s[I, ] <- c(naive.est$Estimate, naive.est$StdError)
  ipw1.s[I, ] <- c(ipw.est1$Estimate, ipw.est1$StdError)
  ipw2.s[I, ] <- c(ipw.est2$Estimate, ipw.est2$StdError)
  ipw3.s[I, ] <- c(ipw.est3$Estimate, ipw.est3$StdError)
  heckman1.s[I, ] <- c(heck.est1$Estimate, heck.est1$StdError)
  heckman2.s[I, ] <- c(heck.est2$Estimate, heck.est2$StdError)
  heckman3.s[I, ] <- c(heck.est3$Estimate, heck.est3$StdError)
  heckman4.s[I, ] <- c(heck.est4$Estimate, heck.est4$StdError)
  heckman5.s[I, ] <- c(heck.est5$Estimate, heck.est5$StdError)
  heckman6.s[I, ] <- c(heck.est6$Estimate, heck.est6$StdError)
  partial1.s[I, ] <- c(partial.est1$Estimate, partial.est1$StdError)
  full.s[I, ] <- c(full.est$Estimate, full.est$StdError)
  partial2.s[I, ] <- c(partial.est2$Estimate, partial.est2$StdError)
  partial3.s[I, ] <- c(partial.est3$Estimate, partial.est3$StdError)
  
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
  diagnostics[I, 11] <- mean(R1[Z1 > 0]) - mean(R1[Z1 < 0])
  diagnostics[I, 12] <- (mean(R1[Z1 > 0]) * (1 - mean(R1[Z1 < 0]))) / (mean(R1[Z1 < 0]) * (1 - mean(R1[Z1 > 0])))
  diagnostics[I, 13] <- mean(R2[Z2 > 0]) - mean(R2[Z2 < 0])
  diagnostics[I, 14] <- (mean(R2[Z2 > 0]) * (1 - mean(R2[Z2 < 0]))) / (mean(R2[Z2 < 0]) * (1 - mean(R2[Z2 > 0])))
  diagnostics[I, 15] <- summary(lm(X1 ~ G1))$f[1]
  diagnostics[I, 16] <- sum(beta.x^2 * 2 * eaf * (1 - eaf)) / var(X1)
  diagnostics[I, 17] <- beta.y^2 * sum(beta.x^2 * 2 * eaf * (1 - eaf)) / var(Y2)
  diagnostics[I, 18] <- summary(lm(X1 ~ G1))$r.sq
  diagnostics[I, 19] <- summary(lm(Y2 ~ G2))$r.sq
  diagnostics[I, 20] <- 1 - glm(R1 ~ Z1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 21] <- 1 - glm(R2 ~ Z2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 22] <- 1 - glm(R1 ~ X1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 23] <- 1 - glm(R2 ~ Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance
  diagnostics[I, 24] <- 1 - glm(R1 ~ X1 + Y1, family = binomial)$deviance / glm(R1 ~ 1, family = binomial)$deviance
  diagnostics[I, 25] <- 1 - glm(R2 ~ X2 + Y2, family = binomial)$deviance / glm(R2 ~ 1, family = binomial)$deviance

  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 20 == 0) save(gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n1, n2, eafs, betas, 
                          diagnostics, bx.oracle.all, sx.oracle.all,
                          bx.naive.all, sx.naive.all, bx.ipw1.all, sx.ipw1.all, bx.ipw2.all, 
                          sx.ipw2.all, bx.ipw3.all, sx.ipw3.all,  
                          bx.heckman1.all, sx.heckman1.all, bx.heckman2.all, sx.heckman2.all, 
                          bx.heckman3.all, sx.heckman3.all, bx.heckman4.all, sx.heckman4.all, 
                          bx.partial1.all, sx.partial1.all, bx.partial2.all, sx.partial2.all, 
                          bx.partial3.all, sx.partial3.all, bx.full.all, sx.full.all, 
                          by.oracle.all, sy.oracle.all, by.naive.all, sy.naive.all, 
                          by.ipw1.all, sy.ipw1.all, by.ipw2.all, sy.ipw2.all, by.ipw3.all, sy.ipw3.all,  
                          by.heckman1.all, sy.heckman1.all, by.heckman2.all, sy.heckman2.all, 
                          by.heckman3.all, sy.heckman3.all, by.heckman4.all, sy.heckman4.all, 
                          by.partial1.all, sy.partial1.all, by.partial2.all, sy.partial2.all, 
                          by.partial3.all, sy.partial3.all, by.full.all, sy.full.all, 
                          oracle.s, naive.s, ipw1.s, ipw2.s, ipw3.s, 
                          heckman1.s, heckman2.s, heckman3.s, heckman4.s, 
                          partial1.s, partial2.s, partial3.s, full.s,
                          bx.heckman5.all, bx.heckman6.all, sx.heckman5.all, 
                          sx.heckman6.all, heckman5.s, heckman6.s, sy.heckman6.all, 
                          by.heckman5.all, by.heckman6.all, sy.heckman5.all, 
                          file = filename)
  
}
