
##########   MR SIMULATION 15   ##########

## Many Instruments - One Sample.

## Load the R functions we will use.
load("IVsel_Functions.RData")

## Load the R package that implements Heckman's method.
library(sampleSelection)
library(ivreg)
library(MendelianRandomization)
library(truncnorm)

## Results will be saved here.
filename <- "MR15_results.RData"

## Set up the data generation.
seed <- 67213383
iter <- 1000
n <- 10000
K <- 10
n.boot1 <- 100
n.boot2 <- 10

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

## Store 2SLS estimates here.
oracle.t <- matrix(0, iter, 4); colnames(oracle.t) <- c("Intercept", "se.int", "Beta", "se.beta")
naive.t <- oracle.t; ipw1.t <- oracle.t; ipw2.t <- oracle.t; ipw3.t <- oracle.t
partial.t <- oracle.t; full.t <- oracle.t; heckman1.t <- oracle.t
heckman2.t <- oracle.t; heckman3.t <- oracle.t; heckman4.t <- oracle.t

## Store the bootstrap samples here.
partial.t.boot <- matrix(0, iter, n.boot2); full.t.boot <- partial.t.boot
heckman1.t.boot <- matrix(0, iter, n.boot1); heckman2.t.boot <- heckman1.t.boot
heckman3.t.boot <- heckman1.t.boot; heckman4.t.boot <- heckman1.t.boot

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
diagnostics <- matrix(0, iter, 23)
colnames(diagnostics) <- c("Strength", "Strength | X", "Strength | X, Y", "Prop Selected", "MD X", "OR X", "MD Y", "OR Y", 
                           "Wmin | X", "Wmax | X", "Wmin | X, Z", "Wmax | X, Z", "MD Z", "OR Z", "F-stat G", "GV% X", "GV% Y",
                           "R2 (X ~ G)", "R2 (Y ~ G)", "McFadden Z", "McFadden X", "McFadden Y", "McFadden XY")

## Start the loop.
for (I in 1:iter) {
  
  ## Seed it.
  set.seed(seed + I)
  
  ## Simulate betas and effect allele frequencies.
  beta.x <- rtruncnorm(K, a = 0.15, mean = 0, sd = 0.05)
  eaf <- runif(K, 0.1, 0.9)
  
  ## Simulate the risk factor, outcome and instruments.
  U <- rnorm(n, 0, 1)
  Z <- rnorm(n, 0, 1)
  G <- matrix(0, n, K)
  for (j in 1:K) G[, j] <- rbinom(n, 2, eaf[j])
  X <- as.vector(G %*% beta.x) + gamma.x * U + rnorm(n, 0, 1)
  Y <- beta.y * X + gamma.y * U + rnorm(n, 0, 1)
  
  ## Simulate the selection coefficient.
  R.probs <- expit(alpha.r + beta.r * X + gamma.r * Z + delta.r * Y)
  R <- rbinom(n, 1, R.probs)
  
  ## Two-stage least squares estimates.
  
  ## Naive and oracle estimates.
  oracle.est <- summary(ivreg(Y ~ X | G))
  Xobs <- X; Xobs[R == 0] <- NA
  Yobs <- Y; Yobs[R == 0] <- NA
  naive.est <- summary(ivreg(Yobs ~ Xobs | G))
  
  ## 2SLS with IPW.
  ipw.est1 <- ipw.ivreg(G, X, Y, R, Z, selection = "G")
  ipw.est2 <- ipw.ivreg(G, X, Y, R, Z, selection = "Z")
  ipw.est3 <- ipw.ivreg(G, X, Y, R, Z, selection = "GZ")
  ## Cannot use either X or Y in the selection model.
  
  ## Heckman model estimates.
  heck.stage1.est1 <- selection(R ~ Z, Xobs ~ G, method = "2step")
  heck.stage1.est2 <- selection(R ~ Z, Xobs ~ G, method = "ml")
  heck.stage1.est3 <- selection(R ~ G + Z, Xobs ~ G, method = "2step")
  heck.stage1.est4 <- selection(R ~ G + Z, Xobs ~ G, method = "ml")
  heck.stage2.est1 <- summary(selection(R ~ Z, Yobs ~ fitted(heck.stage1.est1), method = "2step"))
  heck.stage2.est2 <- summary(selection(R ~ Z, Yobs ~ fitted(heck.stage1.est2), method = "ml"))
  heck.stage2.est3 <- summary(selection(R ~ G + Z, Yobs ~ fitted(heck.stage1.est3), method = "2step"))
  heck.stage2.est4 <- summary(selection(R ~ G + Z, Yobs ~ fitted(heck.stage1.est4), method = "ml"))
  
  ## Bootstrap.
  heck.boot.est1 <- rep(0, n.boot1); heck.boot.est2 <- rep(0, n.boot1)
  heck.boot.est3 <- rep(0, n.boot1); heck.boot.est4 <- rep(0, n.boot1)
  for(i in 1:n.boot1) {
    boot <- sample(n, n, replace = TRUE)
    Gb <- G[boot, ]; Xb <- X[boot]; Rb <- R[boot]; Yb <- Y[boot]; Zb <- Z[boot]
    Xbobs <- Xb; Xbobs[Rb == 0] <- NA; Ybobs <- Yb; Ybobs[Rb == 0] <- NA
    heck.stage1.est1 <- selection(Rb ~ Zb, Xbobs ~ Gb, method = "2step")
    heck.stage1.est2 <- selection(Rb ~ Zb, Xbobs ~ Gb, method = "ml")
    heck.stage1.est3 <- selection(Rb ~ Gb + Zb, Xbobs ~ Gb, method = "2step")
    heck.stage1.est4 <- selection(Rb ~ Gb + Zb, Xbobs ~ Gb, method = "ml")
    heck.boot.est1[i] <- summary(selection(Rb ~ Zb, Ybobs ~ fitted(heck.stage1.est1), method = "2step"))[[2]][4, 1]
    heck.boot.est2[i] <- summary(selection(Rb ~ Zb, Ybobs ~ fitted(heck.stage1.est2), method = "ml"))$estimate[4, 1]
    heck.boot.est3[i] <- summary(selection(Rb ~ Gb + Zb, Ybobs ~ fitted(heck.stage1.est3), method = "2step"))[[2]][K + 4, 1]
    heck.boot.est4[i] <- summary(selection(Rb ~ Gb + Zb, Ybobs ~ fitted(heck.stage1.est4), method = "ml"))$estimate[K + 4, 1]
  }
  
  ## TTW estimates.
  partial.opt1 <- optim(par = rep(0, K + 2), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z)
  partial.opt2 <- optim(par = rep(0, 2 * K + 3), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = X, alpha.hat = partial.opt1$par)
  full.opt <- optim(par = rep(0, 3 * K + 5), full.lik, method = "BFGS", control = list(fnscale = -1), X = G, R = R, Z = Z, Y = X)
  partial.stage1 <- unname( cbind(1, G) %*% partial.opt2$par[1:(K + 1)] )
  full.stage1 <- unname( cbind(1, G) %*% full.opt$par[1:(K + 1)] )
  partial.opt3 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = partial.stage1, R = R, Z = Z)
  partial.est <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = partial.stage1, R = R, Z = Z, Y = Y, alpha.hat = partial.opt3$par, hessian = TRUE)
  partial.sd <- sqrt(diag(solve(- partial.est$hessian))[1:2])
  full.est <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = full.stage1, R = R, Z = Z, Y = Y, hessian = TRUE)
  full.sd <- sqrt(diag(solve(- full.est$hessian))[1:2])
  
  ## Bootstrap.
  ivsel.boot.est1 <- rep(0, n.boot2); ivsel.boot.est2 <- rep(0, n.boot2)
  for(i in 1:n.boot2) {
    boot <- sample(n, n, replace = TRUE)
    Gb <- G[boot, ]; Xb <- X[boot]; Rb <- R[boot]; Yb <- Y[boot]; Zb <- Z[boot]
    partial.opt1 <- optim(par = rep(0, K + 2), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = Gb, R = Rb, Z = Zb)
    partial.opt2 <- optim(par = rep(0, 2 * K + 3), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = Gb, R = Rb, Z = Zb, Y = Xb, alpha.hat = partial.opt1$par)
    full.opt <- optim(par = rep(0, 3 * K + 5), full.lik, method = "BFGS", control = list(fnscale = -1), X = Gb, R = Rb, Z = Zb, Y = Xb)
    partial.stage1 <- unname( cbind(1, Gb) %*% partial.opt2$par[1:(K + 1)] )
    full.stage1 <- unname( cbind(1, Gb) %*% full.opt$par[1:(K + 1)] )
    partial.opt3 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = partial.stage1, R = R, Z = Z)
    ivsel.boot.est1[i] <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = partial.stage1, R = R, Z = Z, Y = Y, alpha.hat = partial.opt3$par)$par[2]
    ivsel.boot.est2[i] <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = full.stage1, R = R, Z = Z, Y = Y)$par[2]
  }
  
  ## Store betas and allele frequencies.
  eafs[I, ] <- eaf
  betas[I, ] <- beta.x
  
  ## Store 2SLS estimates.
  oracle.t[I, ] <- c(oracle.est$coef[1, 1], oracle.est$coef[1, 2], oracle.est$coef[2, 1], oracle.est$coef[2, 2])
  naive.t[I, ] <- c(naive.est$coef[1, 1], naive.est$coef[1, 2], naive.est$coef[2, 1], naive.est$coef[2, 2])
  ipw1.t[I, ] <- c(ipw.est1$est[1, 1], ipw.est1$est[1, 2], ipw.est1$est[2, 1], ipw.est1$est[2, 2])
  ipw2.t[I, ] <- c(ipw.est2$est[1, 1], ipw.est2$est[1, 2], ipw.est2$est[2, 1], ipw.est2$est[2, 2])
  ipw3.t[I, ] <- c(ipw.est3$est[1, 1], ipw.est3$est[1, 2], ipw.est3$est[2, 1], ipw.est3$est[2, 2])
  partial.t[I, ] <- c(partial.est$par[1], partial.sd[1], partial.est$par[2], partial.sd[2])
  full.t[I, ] <- c(full.est$par[1], full.sd[1], full.est$par[2], full.sd[2])
  heckman1.t[I, ] <- c(heck.stage2.est1$estimate[3, 1:2], heck.stage2.est1$estimate[4, 1:2])
  heckman2.t[I, ] <- c(heck.stage2.est2$estimate[3, 1:2], heck.stage2.est2$estimate[4, 1:2])
  heckman3.t[I, ] <- c(heck.stage2.est3$estimate[K + 3, 1:2], heck.stage2.est3$estimate[K + 4, 1:2])
  heckman4.t[I, ] <- c(heck.stage2.est4$estimate[K + 3, 1:2], heck.stage2.est4$estimate[K + 4, 1:2])
  
  ## Store bootstrap samples.
  partial.t.boot[I, ] <- ivsel.boot.est1
  full.t.boot[I, ] <- ivsel.boot.est2
  heckman1.t.boot[I, ] <- heck.boot.est1
  heckman2.t.boot[I, ] <- heck.boot.est2
  heckman3.t.boot[I, ] <- heck.boot.est3
  heckman4.t.boot[I, ] <- heck.boot.est4
  
  ## Store diagnostics.
  diagnostics[I, 9] <- ipw.est1$w.min
  diagnostics[I, 10] <- ipw.est1$w.max
  diagnostics[I, 11] <- ipw.est2$w.min
  diagnostics[I, 12] <- ipw.est2$w.max
  
  
  ## Summary statistics and IVW estimates.
  
  ## Oracle summary statistics.
  bx.oracle <- rep(0, K); sx.oracle <- rep(0, K)
  by.oracle <- rep(0, K); sy.oracle <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(lm(X ~ G[, i]))$coefficients
    bx.oracle[i] <- fit1[2, 1]
    sx.oracle[i] <- fit1[2, 2]
    fit2 <- summary(lm(Y ~ G[, i]))$coefficients
    by.oracle[i] <- fit2[2, 1]
    sy.oracle[i] <- fit2[2, 2]
  }
  
  ## Oracle causal effect estimate.
  oracle.est <- mr_ivw(mr_input(bx = bx.oracle, bxse = sx.oracle, by = by.oracle, byse = sy.oracle))
  
  ## Complete-case summary statistics.
  bx.naive <- rep(0, K); sx.naive <- rep(0, K)
  by.naive <- rep(0, K); sy.naive <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(lm(X[R == 1] ~ G[R == 1, i]))$coefficients
    bx.naive[i] <- fit1[2, 1]
    sx.naive[i] <- fit1[2, 2]
    fit2 <- summary(lm(Y[R == 1] ~ G[R == 1, i]))$coefficients
    by.naive[i] <- fit2[2, 1]
    sy.naive[i] <- fit2[2, 2]
  }
  
  ## Complete-case causal effect estimate.
  naive.est <- mr_ivw(mr_input(bx = bx.naive, bxse = sx.naive, by = by.naive, byse = sy.naive))
  
  ## IPW (G) summary statistics.
  bx.ipw1 <- rep(0, K); sx.ipw1 <- rep(0, K)
  by.ipw1 <- rep(0, K); sy.ipw1 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G[, i], Y = X, R = R, Z = G[, which(!(1:K %in% i))])
    bx.ipw1[i] <- fit1$est[2, 1]
    sx.ipw1[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G[, i], Y = Y, R = R, Z = G[, which(!(1:K %in% i))])
    by.ipw1[i] <- fit2$est[2, 1]
    sy.ipw1[i] <- fit2$est[2, 2]
  }
  
  ## IPW (G) causal effect estimate.
  ipw.est1 <- mr_ivw(mr_input(bx = bx.ipw1, bxse = sx.ipw1, by = by.ipw1, byse = sy.ipw1))
  
  ## IPW (G, Z) summary statistics.
  bx.ipw2 <- rep(0, K); sx.ipw2 <- rep(0, K)
  by.ipw2 <- rep(0, K); sy.ipw2 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G[, i], Y = X, R = R, Z = cbind(Z, G[, which(!(1:K %in% i))]))
    bx.ipw2[i] <- fit1$est[2, 1]
    sx.ipw2[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G[, i], Y = Y, R = R, Z = cbind(Z, G[, which(!(1:K %in% i))]))
    by.ipw2[i] <- fit2$est[2, 1]
    sy.ipw2[i] <- fit2$est[2, 2]
  }
  
  ## IPW (G, Z) causal effect estimate.
  ipw.est2 <- mr_ivw(mr_input(bx = bx.ipw2, bxse = sx.ipw2, by = by.ipw2, byse = sy.ipw2))
  
  ## IPW (Z) summary statistics.
  bx.ipw3 <- rep(0, K); sx.ipw3 <- rep(0, K)
  by.ipw3 <- rep(0, K); sy.ipw3 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- ipw.linear(X = G[, i], Y = X, R = R, Z = Z, drop.x = TRUE)
    bx.ipw3[i] <- fit1$est[2, 1]
    sx.ipw3[i] <- fit1$est[2, 2]
    fit2 <- ipw.linear(X = G[, i], Y = Y, R = R, Z = Z, drop.x = TRUE)
    by.ipw3[i] <- fit2$est[2, 1]
    sy.ipw3[i] <- fit2$est[2, 2]
  }
  
  ## IPW (Z) causal effect estimate.
  ipw.est3 <- mr_ivw(mr_input(bx = bx.ipw3, bxse = sx.ipw3, by = by.ipw3, byse = sy.ipw3))
  
  ## Heckman-2SLS (Z) summary statistics.
  bx.heck1 <- rep(0, K); sx.heck1 <- rep(0, K)
  by.heck1 <- rep(0, K); sy.heck1 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ Z, Xobs ~ G[, i], method = "2step"))[[2]]
    bx.heck1[i] <- fit1[4, 1]
    sx.heck1[i] <- fit1[4, 2]
    fit2 <- summary(selection(R ~ Z, Yobs ~ G[, i], method = "2step"))[[2]]
    by.heck1[i] <- fit2[4, 1]
    sy.heck1[i] <- fit2[4, 2]
  }
  
  ## Heckman-2SLS (Z) causal effect estimate.
  heck.est1 <- mr_ivw(mr_input(bx = bx.heck1, bxse = sx.heck1, by = by.heck1, byse = sy.heck1))
  
  ## Heckman-ML (Z) summary statistics.
  bx.heck2 <- rep(0, K); sx.heck2 <- rep(0, K)
  by.heck2 <- rep(0, K); sy.heck2 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ Z, Xobs ~ G[, i], method = "ml"))$estimate
    bx.heck2[i] <- fit1[4, 1]
    sx.heck2[i] <- fit1[4, 2]
    fit2 <- summary(selection(R ~ Z, Yobs ~ G[, i], method = "ml"))$estimate
    by.heck2[i] <- fit2[4, 1]
    sy.heck2[i] <- fit2[4, 2]
  }
  
  ## Heckman-ML (Z) causal effect estimate.
  heck.est2 <- mr_ivw(mr_input(bx = bx.heck2, bxse = sx.heck2, by = by.heck2, byse = sy.heck2))
  
  ## Heckman-2SLS (G, Z) summary statistics.
  bx.heck3 <- rep(0, K); sx.heck3 <- rep(0, K)
  by.heck3 <- rep(0, K); sy.heck3 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ G + Z, Xobs ~ G[, i], method = "2step"))[[2]]
    bx.heck3[i] <- fit1[K + 4, 1]
    sx.heck3[i] <- fit1[K + 4, 2]
    fit2 <- summary(selection(R ~ G + Z, Yobs ~ G[, i], method = "2step"))[[2]]
    by.heck3[i] <- fit2[K + 4, 1]
    sy.heck3[i] <- fit2[K + 4, 2]
  }
  
  ## Heckman-2SLS (G, Z) causal effect estimate.
  heck.est3 <- mr_ivw(mr_input(bx = bx.heck3, bxse = sx.heck3, by = by.heck3, byse = sy.heck3))
  
  ## Heckman-ML (G, Z) summary statistics.
  bx.heck4 <- rep(0, K); sx.heck4 <- rep(0, K)
  by.heck4 <- rep(0, K); sy.heck4 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ G + Z, Xobs ~ G[, i], method = "ml"))$estimate
    bx.heck4[i] <- fit1[K + 4, 1]
    sx.heck4[i] <- fit1[K + 4, 2]
    fit2 <- summary(selection(R ~ G + Z, Yobs ~ G[, i], method = "ml"))$estimate
    by.heck4[i] <- fit2[K + 4, 1]
    sy.heck4[i] <- fit2[K + 4, 2]
  }
  
  ## Heckman-ML (G, Z) causal effect estimate.
  heck.est4 <- mr_ivw(mr_input(bx = bx.heck4, bxse = sx.heck4, by = by.heck4, byse = sy.heck4))
  
  ## Heckman-2SLS (G, Z) summary statistics.
  bx.heck5 <- rep(0, K); sx.heck5 <- rep(0, K)
  by.heck5 <- rep(0, K); sy.heck5 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ G[, i] + Z, Xobs ~ G[, i], method = "2step"))[[2]]
    bx.heck5[i] <- fit1[5, 1]
    sx.heck5[i] <- fit1[5, 2]
    fit2 <- summary(selection(R ~ G[, i] + Z, Yobs ~ G[, i], method = "2step"))[[2]]
    by.heck5[i] <- fit2[5, 1]
    sy.heck5[i] <- fit2[5, 2]
  }
  
  ## Heckman-2SLS (G, Z) causal effect estimate.
  heck.est5 <- mr_ivw(mr_input(bx = bx.heck5, bxse = sx.heck5, by = by.heck5, byse = sy.heck5))
  
  ## Heckman-ML (G, Z) summary statistics.
  bx.heck6 <- rep(0, K); sx.heck6 <- rep(0, K)
  by.heck6 <- rep(0, K); sy.heck6 <- rep(0, K)
  for (i in 1:K) {
    fit1 <- summary(selection(R ~ G[, i] + Z, Xobs ~ G[, i], method = "ml"))$estimate
    bx.heck6[i] <- fit1[5, 1]
    sx.heck6[i] <- fit1[5, 2]
    fit2 <- summary(selection(R ~ G[, i] + Z, Yobs ~ G[, i], method = "ml"))$estimate
    by.heck6[i] <- fit2[5, 1]
    sy.heck6[i] <- fit2[5, 2]
  }
  
  ## Heckman-ML (G, Z) causal effect estimate.
  heck.est6 <- mr_ivw(mr_input(bx = bx.heck6, bxse = sx.heck6, by = by.heck6, byse = sy.heck6))
  
  ## TTW-partial (default) summary statistics.
  bx.partial1 <- rep(0, K); sx.partial1 <- rep(0, K)
  by.partial1 <- rep(0, K); sy.partial1 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z)
    partial.opt2 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = X, alpha.hat = partial.opt1$par, hessian = TRUE)
    bx.partial1[i] <- partial.opt2$par[2]
    sx.partial1[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt3 <- optim(par = rep(0, 5), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = Y, alpha.hat = partial.opt1$par, hessian = TRUE)
    by.partial1[i] <- partial.opt3$par[2]
    sy.partial1[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (default) causal effect estimate.
  partial.est1 <- mr_ivw(mr_input(bx = bx.partial1, bxse = sx.partial1, by = by.partial1, byse = sy.partial1))
  
  ## TTW-full (default) summary statistics.
  bx.full <- rep(0, K); sx.full <- rep(0, K)
  by.full <- rep(0, K); sy.full <- rep(0, K)
  for (i in 1:K) {
    full.opt1 <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = X, hessian = TRUE)
    bx.full[i] <- full.opt1$par[2]
    sx.full[i] <- sqrt(diag(solve(- full.opt1$hessian))[2])
    full.opt2 <- optim(par = rep(0, 8), full.lik, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = Y, hessian = TRUE)
    by.full[i] <- full.opt2$par[2]
    sy.full[i] <- sqrt(diag(solve(- full.opt2$hessian))[2])
  }
  
  ## TTW-full (default) causal effect estimate.
  full.est <- mr_ivw(mr_input(bx = bx.full, bxse = sx.full, by = by.full, byse = sy.full))
  
  ## TTW-partial (G-prop) summary statistics.
  bx.partial2 <- rep(0, K); sx.partial2 <- rep(0, K)
  by.partial2 <- rep(0, K); sy.partial2 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1 <- optim(par = rep(0, 12), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = cbind(Z, G[, which(!(1:K %in% i))]))
    partial.opt2 <- optim(par = rep(0, 5), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = cbind(Z, G[, which(!(1:K %in% i))]), Y = X, alpha.hat = partial.opt1$par, hessian = TRUE)
    bx.partial2[i] <- partial.opt2$par[2]
    sx.partial2[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt3 <- optim(par = rep(0, 5), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = cbind(Z, G[, which(!(1:K %in% i))]), Y = Y, alpha.hat = partial.opt1$par, hessian = TRUE)
    by.partial2[i] <- partial.opt3$par[2]
    sy.partial2[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (G-prop) causal effect estimate.
  partial.est2 <- mr_ivw(mr_input(bx = bx.partial2, bxse = sx.partial2, by = by.partial2, byse = sy.partial2))
  
  ## TTW-partial (G-sel) summary statistics.
  bx.partial3 <- rep(0, K); sx.partial3 <- rep(0, K)
  by.partial3 <- rep(0, K); sy.partial3 <- rep(0, K)
  for (i in 1:K) {
    partial.opt1 <- optim(par = rep(0, 3), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z)
    partial.opt2 <- optim(par = rep(0, 14), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = X, add.sel = G[, which(!(1:K %in% i))], alpha.hat = partial.opt1$par, hessian = TRUE)
    bx.partial3[i] <- partial.opt2$par[2]
    sx.partial3[i] <- sqrt(diag(solve(- partial.opt2$hessian))[2])
    partial.opt3 <- optim(par = rep(0, 14), partial.lik2.v2, method = "BFGS", control = list(fnscale = -1), X = G[, i], R = R, Z = Z, Y = Y, add.sel = G[, which(!(1:K %in% i))], alpha.hat = partial.opt1$par, hessian = TRUE)
    by.partial3[i] <- partial.opt3$par[2]
    sy.partial3[i] <- sqrt(diag(solve(- partial.opt3$hessian))[2])
  }
  
  ## TTW-partial (G-sel) causal effect estimate.
  partial.est3 <- mr_ivw(mr_input(bx = bx.partial3, bxse = sx.partial3, by = by.partial3, byse = sy.partial3))
  
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
  diagnostics[I, 1] <- glm(R ~ 1, family = binomial)$deviance - glm(R ~ Z, family = binomial)$deviance
  diagnostics[I, 2] <- glm(R ~ X, family = binomial)$deviance - glm(R ~ X + Z, family = binomial)$deviance
  diagnostics[I, 3] <- glm(R ~ X + Y, family = binomial)$deviance - glm(R ~ X + Y + Z, family = binomial)$deviance
  diagnostics[I, 4] <- mean(R)
  diagnostics[I, 5] <- mean(R[X > 0]) - mean(R[X < 0])
  diagnostics[I, 6] <- (mean(R[X > 0]) * (1 - mean(R[X < 0]))) / (mean(R[X < 0]) * (1 - mean(R[X > 0])))
  diagnostics[I, 7] <- mean(R[Y > 1]) - mean(R[Y < 1])
  diagnostics[I, 8] <- (mean(R[Y > 1]) * (1 - mean(R[Y < 1]))) / (mean(R[Y < 1]) * (1 - mean(R[Y > 1])))
  diagnostics[I, 13] <- mean(R[Z > 0]) - mean(R[Z < 0])
  diagnostics[I, 14] <- (mean(R[Z > 0]) * (1 - mean(R[Z < 0]))) / (mean(R[Z < 0]) * (1 - mean(R[Z > 0])))
  diagnostics[I, 15] <- summary(lm(X ~ G))$f[1]
  diagnostics[I, 16] <- sum(beta.x^2 * 2 * eaf * (1 - eaf)) / var(X)
  diagnostics[I, 17] <- beta.y^2 * sum(beta.x^2 * 2 * eaf * (1 - eaf)) / var(Y)
  diagnostics[I, 18] <- summary(lm(X ~ G))$r.sq
  diagnostics[I, 19] <- summary(lm(Y ~ G))$r.sq
  diagnostics[I, 20] <- 1 - glm(R ~ Z, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 21] <- 1 - glm(R ~ X, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 22] <- 1 - glm(R ~ Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  diagnostics[I, 23] <- 1 - glm(R ~ X + Y, family = binomial)$deviance / glm(R ~ 1, family = binomial)$deviance
  
  ## Print progress status.
  if (I %% 20 == 0) print(paste(I, " done."))
  
  ## Save results.
  if (I %% 20 == 0) save(gamma.x, beta.y, gamma.y, alpha.r, beta.r, gamma.r, delta.r, n, eafs, betas, 
                          oracle.t, naive.t, ipw1.t, ipw2.t, ipw3.t, diagnostics,
                          heckman1.t, heckman2.t, heckman3.t, heckman4.t, partial.t, full.t, 
                          heckman1.t.boot, heckman2.t.boot, heckman3.t.boot, heckman4.t.boot, 
                          partial.t.boot, full.t.boot, bx.oracle.all, sx.oracle.all,
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
