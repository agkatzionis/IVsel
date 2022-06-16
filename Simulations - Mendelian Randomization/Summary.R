
##########   SUMMARY   ##########

## Here we assess the performance of IVsel methods in MR  
## simulations and create the Tables included in the paper.

## Set working directory.
setwd( "MY WORKING DIRECTORY" )

## Load R packages and relevant functions.
load("IVsel_Functions.RData")
library(sampleSelection)
library(ivreg)
library(MendelianRandomization)
library(truncnorm)
library(xtable)

##################################################

## Simulation results are supplied in separate files. 
## Here, we combine those files to create Tables 3 and 4,
## as well as Supplementary Tables 5 and 6 of the paper.


## In each Table we report causal effect estimates, the
## standard deviation of estimates, (second-order) causal 
## standard errors, coverage/Type I error rates and  (for 
## theta != 0) the power to reject the causal null hypothesis.

##################################################

##########   TABLE 3   ##########

## This Table contains results of implementing the IVsel methods 
## in MR simulations with a single genetic instrument for inference.
## Results are contained in the files "MR1-12_results.RData".

## Create the table.
Table3 <- matrix(0, 30, 9)
colnames(Table3) <- c("Causal", " Emp SE", "Causal SE", "Coverage", "Power", "Causal", "Emp SE", "Causal SE", "Type I")
rownames(Table3) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), times = 6)

## The 12 scenarios are for X, Y, X&Y affecting selection respectively,
## for either a positive (theta = 0.2) or a null (theta = 0) causal 
## effect, first in one-sample and then in two-sample MR.

## Scenario 1: One-sample MR, theta = 0.2, Y --> R.
load("MR1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[1:5, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[1:5, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[1:5, 3] <- stderr
Table3[1:5, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[1:5, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2: One-sample MR, theta = 0.2, X --> R.
load("MR2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / ipw1[, 3]^2 + gy[, 3]^2 * ipw1[, 4]^2 / ipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[6:10, 1] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / ipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[6:10, 2] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / ipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[6:10, 3] <- stderr
Table3[6:10, 4] <- c(mean( ((gy[, 3] / naive[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[6:10, 5] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3: One-sample MR, theta = 0.2, X,Y --> R.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 3] <- stderr
Table3[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4: One-sample MR, theta = 0, Y --> R.
load("MR4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[1:5, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[1:5, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[1:5, 8] <- stderr
Table3[1:5, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5: One-sample MR, theta = 0, X --> R.
load("MR5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / ipw1[, 3]^2 + gy[, 3]^2 * ipw1[, 4]^2 / ipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[6:10, 6] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / ipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[6:10, 7] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / ipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[6:10, 8] <- stderr
Table3[6:10, 9] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6: One-sample MR, theta = 0, X,Y --> R.
load("MR6_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[11:15, 6] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 7] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 8] <- stderr
Table3[11:15, 9] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7: Two-sample MR, theta = 0.2, Y --> R.
load("MR7_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[16:20, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[16:20, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[16:20, 3] <- stderr
Table3[16:20, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[16:20, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 8: Two-sample MR, theta = 0.2, X --> R.
load("MR8_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / ipw1[, 3]^2 + gy[, 3]^2 * ipw1[, 4]^2 / ipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[21:25, 1] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / ipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[21:25, 2] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / ipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[21:25, 3] <- stderr
Table3[21:25, 4] <- c(mean( ((gy[, 3] / naive[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[21:25, 5] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 9: Two-sample MR, theta = 0.2, X,Y --> R.
load("MR9_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 3] <- stderr
Table3[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 10: Two-sample MR, theta = 0, Y --> R.
load("MR10_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[16:20, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[16:20, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[16:20, 8] <- stderr
Table3[16:20, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 11: Two-sample MR, theta = 0, X --> R.
load("MR11_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / ipw1[, 3]^2 + gy[, 3]^2 * ipw1[, 4]^2 / ipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[21:25, 6] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / ipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[21:25, 7] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / ipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[21:25, 8] <- stderr
Table3[21:25, 9] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / ipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / ipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 12: Two-sample MR, theta = 0, XY --> R.
load("MR12_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[26:30, 6] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 7] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 8] <- stderr
Table3[26:30, 9] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   TABLE 4   ##########

## This Table contains results of implementing the IVsel methods 
## in MR simulations with multiple instruments for inference.
## Results are contained in the files "MR13-24_results.RData".

## Create the table.
Table4 <- matrix(0, 30, 9)
colnames(Table4) <- c("Causal", " Emp SE", "Causal SE", "Coverage", "Power", "Causal", "Emp SE", "Causal SE", "Type I")
rownames(Table4) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), times = 6)

## The 12 scenarios are for X, Y, X&Y affecting selection respectively,
## for either a positive (theta = 0.2) or a null (theta = 0) causal 
## effect, first in one-sample and then in two-sample MR.

## Scenario 13: One-sample MR, theta = 0.2, Y --> R.
load("MR13_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[1:5, 1] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[1:5, 2] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[1:5, 3] <- stderr
Table4[1:5, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (ipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (ipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[1:5, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 14: One-sample MR, theta = 0.2, X --> R.
load("MR14_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[6:10, 1] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[6:10, 2] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[6:10, 3] <- stderr
Table4[6:10, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (ipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (ipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[6:10, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 15: One-sample MR, theta = 0.2, X,Y --> R.
load("MR15_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[11:15, 1] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[11:15, 2] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[11:15, 3] <- stderr
Table4[11:15, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (ipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (ipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[11:15, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 16: One-sample MR, theta = 0, Y --> R.
load("MR16_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[1:5, 6] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[1:5, 7] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[1:5, 8] <- stderr
Table4[1:5, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 17: One-sample MR, theta = 0, X --> R.
load("MR17_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[6:10, 6] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[6:10, 7] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[6:10, 8] <- stderr
Table4[6:10, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 18: One-sample MR, theta = 0, X,Y --> R.
load("MR18_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(ipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[11:15, 6] <- c(mean(naive.t[, 3]), mean(ipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[11:15, 7] <- c(sd(naive.t[, 3]), sd(ipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[11:15, 8] <- stderr
Table4[11:15, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.t[, 3] - 1.96 * stderr[2]) * (ipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 19: Two-sample MR, theta = 0.2, Y --> R.
load("MR19_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[16:20, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[16:20, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[16:20, 3] <- stderr
Table4[16:20, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[16:20, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 20: Two-sample MR, theta = 0.2, X --> R.
load("MR20_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[21:25, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[21:25, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[21:25, 3] <- stderr
Table4[21:25, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[21:25, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 21: Two-sample MR, theta = 0.2, X,Y --> R.
load("MR21_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[26:30, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[26:30, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[26:30, 3] <- stderr
Table4[26:30, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[26:30, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 22: Two-sample MR, theta = 0, Y --> R.
load("MR22_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[16:20, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[16:20, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[16:20, 8] <- stderr
Table4[16:20, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 23: Two-sample MR, theta = 0, X --> R.
load("MR23_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[21:25, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[21:25, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[21:25, 8] <- stderr
Table4[21:25, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 24: Two-sample MR, theta = 0, X,Y --> R.
load("MR24_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[26:30, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[26:30, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[26:30, 8] <- stderr
Table4[26:30, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

##################################################

##########   SUPPLEMENTARY TABLE 5   ##########

## This Table contains results of implementing the summary-statistics
## version of the IVsel methods in one-sample MR simulations 
## with multiple instruments for inference. Results are contained 
## in the files "MR13-18_results.RData".

## Create the table.
SuppTable5 <- matrix(0, 15, 9)
colnames(SuppTable5) <- c("Causal", " Emp SE", "Causal SE", "Coverage", "Power", "Causal", "Emp SE", "Causal SE", "Type I")
rownames(SuppTable5) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), times = 3)

## The 6 scenarios are for X, Y, X&Y affecting selection respectively,
## for either a positive (theta = 0.2) or a null (theta = 0) effect.

## Scenario 13: One-sample MR, theta = 0.2, Y --> R.
load("MR13_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[1:5, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[1:5, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[1:5, 3] <- stderr
SuppTable5[1:5, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[1:5, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 14: One-sample MR, theta = 0.2, X --> R.
load("MR14_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[6:10, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[6:10, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[6:10, 3] <- stderr
SuppTable5[6:10, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[6:10, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 15: One-sample MR, theta = 0.2, X,Y --> R.
load("MR15_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[11:15, 1] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[11:15, 2] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[11:15, 3] <- stderr
SuppTable5[11:15, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (ipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (ipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[11:15, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 16: One-sample MR, theta = 0, Y --> R.
load("MR16_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[1:5, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[1:5, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[1:5, 8] <- stderr
SuppTable5[1:5, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 17: One-sample MR, theta = 0, X --> R.
load("MR17_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[6:10, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[6:10, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[6:10, 8] <- stderr
SuppTable5[6:10, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 18: One-sample MR, theta = 0, X,Y --> R.
load("MR18_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.ipw1.all, sx.ipw1.all, by.ipw1.all, sy.ipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[11:15, 6] <- c(mean(naive.s[, 1]), mean(ipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[11:15, 7] <- c(sd(naive.s[, 1]), sd(ipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[11:15, 8] <- stderr
SuppTable5[11:15, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (ipw1.s[, 1] - 1.96 * stderr[2]) * (ipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

##################################################

##########   SUPPLEMENTARY TABLE 6   ##########

## This Table contains results of implementing the IVsel methods
## in MR simulations with a single genetic instrument for inference, 
## when the first-stage model is potentially misspecified. Two
## of these simulations are copied from Scenarios 1 and 7 in Table 3; 
## the rest are contained in the files "MR25-28_results.RData".

## Create the table.
SuppTable6 <- matrix(0, 15, 10)
colnames(SuppTable6) <- c("Causal", " Emp SE", "Causal SE", "Coverage", "Power", "Causal", "Emp SE", "Causal SE", "Coverage", "Power")
rownames(SuppTable6) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), times = 3)

## The 6 scenarios represent a linear or non-linear (misspecified) 
## first-stage model and either one-sample MR, or two-sample MR with 
## samples from the same population, or two-sample MR with samples from
## different populations. All simulations are for theta = 0.2, Y --> R.

## Scenario 1: One-sample MR, linear 1st-stage.
load("MR1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[1:5, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 3] <- stderr
SuppTable6[1:5, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[1:5, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 25: One-sample MR, non-linear 1st-stage.
load("MR25_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[1:5, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 8] <- stderr
SuppTable6[1:5, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[1:5, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )


## Scenario 7: Two-sample MR, same population, linear 1st-stage.
load("MR7_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[6:10, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                       mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                       sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 3] <- stderr
SuppTable6[6:10, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[6:10, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 26: Two-sample MR, same population, non-linear 1st-stage.
load("MR26_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[6:10, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                       mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                       sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 8] <- stderr
SuppTable6[6:10, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[6:10, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 27: Two-sample MR, different populations, linear 1st-stage.
load("MR27_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[11:15, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                        mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                        sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 3] <- stderr
SuppTable6[11:15, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[11:15, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 28: Two-sample MR, different populations, non-linear 1st-stage.
load("MR28_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(ipw1[, 4]^2 / gx[, 3]^2 + ipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[11:15, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(ipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                        mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(ipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                        sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 8] <- stderr
SuppTable6[11:15, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[11:15, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((ipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((ipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

## Print all tables in LaTeX format.
xtable(Table3, digits = 3)
xtable(Table4, digits = 3)
xtable(SuppTable5, digits = 3)
xtable(SuppTable6, digits = 3)

##################################################
