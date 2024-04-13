
##########   SUMMARY   ##########

## Here we assess the performance of IVsel methods in MR  
## simulations and create the Tables included in the paper.

## Set working directory.
setwd( "MY WORKING DIRECTORY" )

## Load R packages and relevant functions.
load("IVsel_Functions.RData")
library(sampleSelection)
library(survey)
library(ivreg)
library(MendelianRandomization)
library(truncnorm)
library(xtable)

##################################################

## Simulation results are supplied in separate files. 
## Here, we combine those files to create Tables 3 and 4,
## as well as Supplementary Tables 5 and 6 and 
## Supplementary Figures 5 and 6 of the paper.


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
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[1:5, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[1:5, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[1:5, 3] <- stderr
Table3[1:5, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[1:5, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2: One-sample MR, theta = 0.2, X --> R.
load("MR2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / svyipw1[, 3]^2 + gy[, 3]^2 * svyipw1[, 4]^2 / svyipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[6:10, 1] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / svyipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[6:10, 2] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / svyipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[6:10, 3] <- stderr
Table3[6:10, 4] <- c(mean( ((gy[, 3] / naive[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[6:10, 5] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3: One-sample MR, theta = 0.2, X,Y --> R.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(svyipw1.y[, 3] / svyipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 3] <- stderr
Table3[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4: One-sample MR, theta = 0, Y --> R.
load("MR4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[1:5, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[1:5, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[1:5, 8] <- stderr
Table3[1:5, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5: One-sample MR, theta = 0, X --> R.
load("MR5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / svyipw1[, 3]^2 + gy[, 3]^2 * svyipw1[, 4]^2 / svyipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[6:10, 6] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / svyipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[6:10, 7] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / svyipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[6:10, 8] <- stderr
Table3[6:10, 9] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6: One-sample MR, theta = 0, X,Y --> R.
load("MR6_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[11:15, 6] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(svyipw1.y[, 3] / svyipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 7] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[11:15, 8] <- stderr
Table3[11:15, 9] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7: Two-sample MR, theta = 0.2, Y --> R.
load("MR7_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[16:20, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[16:20, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[16:20, 3] <- stderr
Table3[16:20, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[16:20, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 8: Two-sample MR, theta = 0.2, X --> R.
load("MR8_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / svyipw1[, 3]^2 + gy[, 3]^2 * svyipw1[, 4]^2 / svyipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[21:25, 1] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / svyipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[21:25, 2] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / svyipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[21:25, 3] <- stderr
Table3[21:25, 4] <- c(mean( ((gy[, 3] / naive[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[21:25, 5] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 9: Two-sample MR, theta = 0.2, X,Y --> R.
load("MR9_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(svyipw1.y[, 3] / svyipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 3] <- stderr
Table3[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
Table3[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 10: Two-sample MR, theta = 0, Y --> R.
load("MR10_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
Table3[16:20, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
Table3[16:20, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
Table3[16:20, 8] <- stderr
Table3[16:20, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 11: Two-sample MR, theta = 0, X --> R.
load("MR11_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(gy[, 4]^2 / naive[, 3]^2 + gy[, 3]^2 * naive[, 4]^2 / naive[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / svyipw1[, 3]^2 + gy[, 3]^2 * svyipw1[, 4]^2 / svyipw1[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / heckman3[, 3]^2 + gy[, 3]^2 * heckman3[, 4]^2 / heckman3[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / partial[, 3]^2 + gy[, 3]^2 * partial[, 4]^2 / partial[, 3]^4)), 
            mean(sqrt(gy[, 4]^2 / oracle[, 3]^2 + gy[, 3]^2 * oracle[, 4]^2 / oracle[, 3]^4)) )

## Store relevant values.
Table3[21:25, 6] <- c(mean(gy[, 3] / naive[, 3]), mean(gy[, 3] / svyipw1[, 3]), mean(gy[, 3] / heckman3[, 3]), 
                       mean(gy[, 3] / partial[, 3]), mean(gy[, 3] / oracle[, 3]) )
Table3[21:25, 7] <- c(sd(gy[, 3] / naive[, 3]), sd(gy[, 3] / svyipw1[, 3]), sd(gy[, 3] / heckman3[, 3]), 
                       sd(gy[, 3] / partial[, 3]), sd(gy[, 3] / oracle[, 3]) )
Table3[21:25, 8] <- stderr
Table3[21:25, 9] <- c(mean( ((gy[, 3] / naive[, 3]) - 1.96 * (stderr[1])) * ((gy[, 3] / naive[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((gy[, 3] / svyipw1[, 3]) - 1.96 * (stderr[2])) * ((gy[, 3] / svyipw1[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((gy[, 3] / heckman3[, 3]) - 1.96 * (stderr[3])) * ((gy[, 3] / heckman3[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((gy[, 3] / partial[, 3]) - 1.96 * (stderr[4])) * ((gy[, 3] / partial[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((gy[, 3] / oracle[, 3]) - 1.96 * (stderr[5])) * ((gy[, 3] / oracle[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 12: Two-sample MR, theta = 0, XY --> R.
load("MR12_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
Table3[26:30, 6] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(svyipw1.y[, 3] / svyipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 7] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
Table3[26:30, 8] <- stderr
Table3[26:30, 9] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
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
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[1:5, 1] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[1:5, 2] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[1:5, 3] <- stderr
Table4[1:5, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (svyipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (svyipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[1:5, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 14: One-sample MR, theta = 0.2, X --> R.
load("MR14_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[6:10, 1] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[6:10, 2] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[6:10, 3] <- stderr
Table4[6:10, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (svyipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (svyipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[6:10, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 15: One-sample MR, theta = 0.2, X,Y --> R.
load("MR15_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[11:15, 1] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[11:15, 2] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[11:15, 3] <- stderr
Table4[11:15, 4] <- c(mean( (naive.t[, 3] - 0.2 - 1.96 * stderr[1]) * (naive.t[, 3] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (svyipw1.t[, 3] - 0.2 - 1.96 * stderr[2]) * (svyipw1.t[, 3] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman3.t[, 3] - 0.2 - 1.96 * stderr[3]) * (heckman3.t[, 3] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial.t[, 3] - 0.2 - 1.96 * stderr[4]) * (partial.t[, 3] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.t[, 3] - 0.2 - 1.96 * stderr[5]) * (oracle.t[, 3] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[11:15, 5] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 16: One-sample MR, theta = 0, Y --> R.
load("MR16_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[1:5, 6] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[1:5, 7] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[1:5, 8] <- stderr
Table4[1:5, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 17: One-sample MR, theta = 0, X --> R.
load("MR17_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[6:10, 6] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[6:10, 7] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[6:10, 8] <- stderr
Table4[6:10, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 18: One-sample MR, theta = 0, X,Y --> R.
load("MR18_results.RData")

## Compute standard errors (here no second order but bootstrap).
stderr <- c(mean(naive.t[, 4]), mean(svyipw1.t[, 4]), mean(apply(heckman3.t.boot, 1, sd)), 
            mean(apply(partial.t.boot, 1, sd)), mean(oracle.t[, 4]) )

## Store relevant values.
Table4[11:15, 6] <- c(mean(naive.t[, 3]), mean(svyipw1.t[, 3]), mean(heckman3.t[, 3]), mean(partial.t[, 3]), mean(oracle.t[, 3]) )
Table4[11:15, 7] <- c(sd(naive.t[, 3]), sd(svyipw1.t[, 3]), sd(heckman3.t[, 3]), sd(partial.t[, 3]), sd(oracle.t[, 3]) )
Table4[11:15, 8] <- stderr
Table4[11:15, 9] <- c(mean( (naive.t[, 3] - 1.96 * stderr[1]) * (naive.t[, 3] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.t[, 3] - 1.96 * stderr[2]) * (svyipw1.t[, 3] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman3.t[, 3] - 1.96 * stderr[3]) * (heckman3.t[, 3] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial.t[, 3] - 1.96 * stderr[4]) * (partial.t[, 3] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.t[, 3] - 1.96 * stderr[5]) * (oracle.t[, 3] + 1.96 * stderr[5]) > 0 ) )

## Scenario 19: Two-sample MR, theta = 0.2, Y --> R.
load("MR19_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[16:20, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[16:20, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[16:20, 3] <- stderr
Table4[16:20, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[16:20, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 20: Two-sample MR, theta = 0.2, X --> R.
load("MR20_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[21:25, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[21:25, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[21:25, 3] <- stderr
Table4[21:25, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[21:25, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 21: Two-sample MR, theta = 0.2, X,Y --> R.
load("MR21_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[26:30, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[26:30, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[26:30, 3] <- stderr
Table4[26:30, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
Table4[26:30, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 22: Two-sample MR, theta = 0, Y --> R.
load("MR22_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[16:20, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[16:20, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[16:20, 8] <- stderr
Table4[16:20, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 23: Two-sample MR, theta = 0, X --> R.
load("MR23_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[21:25, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[21:25, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[21:25, 8] <- stderr
Table4[21:25, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 24: Two-sample MR, theta = 0, X,Y --> R.
load("MR24_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
Table4[26:30, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
Table4[26:30, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
Table4[26:30, 8] <- stderr
Table4[26:30, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
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
            mean( apply(cbind(bx.all, sx.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[1:5, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[1:5, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[1:5, 3] <- stderr
SuppTable5[1:5, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                      mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                      mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                      mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                      mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[1:5, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 14: One-sample MR, theta = 0.2, X --> R.
load("MR14_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[6:10, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[6:10, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[6:10, 3] <- stderr
SuppTable5[6:10, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                       mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                       mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                       mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                       mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[6:10, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 15: One-sample MR, theta = 0.2, X,Y --> R.
load("MR15_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[11:15, 1] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[11:15, 2] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[11:15, 3] <- stderr
SuppTable5[11:15, 4] <- c(mean( (naive.s[, 1] - 0.2 - 1.96 * stderr[1]) * (naive.s[, 1] - 0.2 + 1.96 * stderr[1]) < 0 ), 
                        mean( (svyipw1.s[, 1] - 0.2 - 1.96 * stderr[2]) * (svyipw1.s[, 1] - 0.2 + 1.96 * stderr[2]) < 0 ), 
                        mean( (heckman5.s[, 1] - 0.2 - 1.96 * stderr[3]) * (heckman5.s[, 1] - 0.2 + 1.96 * stderr[3]) < 0 ), 
                        mean( (partial1.s[, 1] - 0.2 - 1.96 * stderr[4]) * (partial1.s[, 1] - 0.2 + 1.96 * stderr[4]) < 0 ), 
                        mean( (oracle.s[, 1] - 0.2 - 1.96 * stderr[5]) * (oracle.s[, 1] - 0.2 + 1.96 * stderr[5]) < 0 ) )
SuppTable5[11:15, 5] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                        mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                        mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                        mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 16: One-sample MR, theta = 0, Y --> R.
load("MR16_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.all, sx.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.all, sx.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[1:5, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[1:5, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[1:5, 8] <- stderr
SuppTable5[1:5, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                      mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                      mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                      mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                      mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 17: One-sample MR, theta = 0, X --> R.
load("MR17_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.all, sy.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[6:10, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[6:10, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[6:10, 8] <- stderr
SuppTable5[6:10, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                       mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
                       mean( (heckman5.s[, 1] - 1.96 * stderr[3]) * (heckman5.s[, 1] + 1.96 * stderr[3]) > 0 ), 
                       mean( (partial1.s[, 1] - 1.96 * stderr[4]) * (partial1.s[, 1] + 1.96 * stderr[4]) > 0 ), 
                       mean( (oracle.s[, 1] - 1.96 * stderr[5]) * (oracle.s[, 1] + 1.96 * stderr[5]) > 0 ) )

## Scenario 18: One-sample MR, theta = 0, X,Y --> R.
load("MR18_results.RData")

## Compute standard errors (using mr_ivw).
stderr <- c(mean( apply(cbind(bx.naive.all, sx.naive.all, by.naive.all, sy.naive.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.svyipw1.all, sx.svyipw1.all, by.svyipw1.all, sy.svyipw1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.heckman5.all, sx.heckman5.all, by.heckman5.all, sy.heckman5.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.partial1.all, sx.partial1.all, by.partial1.all, sy.partial1.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ),
            mean( apply(cbind(bx.oracle.all, sx.oracle.all, by.oracle.all, sy.oracle.all), 1, function (x) mr_ivw(mr_input(bx = x[1:10], bxse = x[11:20], by = x[21:30], byse = x[31:40]), weights = "delta")$StdError) ) )

## Store relevant values.
SuppTable5[11:15, 6] <- c(mean(naive.s[, 1]), mean(svyipw1.s[, 1]), mean(heckman5.s[, 1]), mean(partial1.s[, 1]), mean(oracle.s[, 1]) )
SuppTable5[11:15, 7] <- c(sd(naive.s[, 1]), sd(svyipw1.s[, 1]), sd(heckman5.s[, 1]), sd(partial1.s[, 1]), sd(oracle.s[, 1]) )
SuppTable5[11:15, 8] <- stderr
SuppTable5[11:15, 9] <- c(mean( (naive.s[, 1] - 1.96 * stderr[1]) * (naive.s[, 1] + 1.96 * stderr[1]) > 0 ), 
                        mean( (svyipw1.s[, 1] - 1.96 * stderr[2]) * (svyipw1.s[, 1] + 1.96 * stderr[2]) > 0 ), 
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
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[1:5, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 3] <- stderr
SuppTable6[1:5, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[1:5, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 25: One-sample MR, non-linear 1st-stage.
load("MR25_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[1:5, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                      mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                      sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[1:5, 8] <- stderr
SuppTable6[1:5, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[1:5, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )


## Scenario 7: Two-sample MR, same population, linear 1st-stage.
load("MR7_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[6:10, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                       mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                       sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 3] <- stderr
SuppTable6[6:10, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[6:10, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 26: Two-sample MR, same population, non-linear 1st-stage.
load("MR26_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[6:10, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                       mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                       sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[6:10, 8] <- stderr
SuppTable6[6:10, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[6:10, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 27: Two-sample MR, different populations, linear 1st-stage.
load("MR27_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[11:15, 1] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                        mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 2] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                        sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 3] <- stderr
SuppTable6[11:15, 4] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[11:15, 5] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 28: Two-sample MR, different populations, non-linear 1st-stage.
load("MR28_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive[, 4]^2 / gx[, 3]^2 + naive[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(svyipw1[, 4]^2 / gx[, 3]^2 + svyipw1[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(heckman3[, 4]^2 / gx[, 3]^2 + heckman3[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)),
            mean(sqrt(partial[, 4]^2 / gx[, 3]^2 + partial[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)), 
            mean(sqrt(oracle[, 4]^2 / gx[, 3]^2 + oracle[, 3]^2 * gx[, 4]^2 / gx[, 3]^4)) )

## Store relevant values.
SuppTable6[11:15, 6] <- c(mean(naive[, 3] / gx[, 3]), mean(svyipw1[, 3] / gx[, 3]), mean(heckman3[, 3] / gx[, 3]),
                        mean(partial[, 3] / gx[, 3]), mean(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 7] <- c(sd(naive[, 3] / gx[, 3]), sd(svyipw1[, 3] / gx[, 3]), sd(heckman3[, 3] / gx[, 3]),
                        sd(partial[, 3] / gx[, 3]), sd(oracle[, 3] / gx[, 3]) )
SuppTable6[11:15, 8] <- stderr
SuppTable6[11:15, 9] <- c(mean( ((naive[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle[, 3] / gx[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle[, 3] / gx[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
SuppTable6[11:15, 10] <- c(mean( ((naive[, 3] / gx[, 3]) - 1.96 * (stderr[1])) * ((naive[, 3] / gx[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1[, 3] / gx[, 3]) - 1.96 * (stderr[2])) * ((svyipw1[, 3] / gx[, 3]) + 1.96 * (stderr[2])) > 0 ), 
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

## Now create the two Figures. Due to high uncertainty
## for the TTW method, we have outliers that mildly 
## imbalance the plots. For better visualization, 
## we plot median values instead of means.

##########   SUPPLEMENTARY PLOT 5   ##########

## To make it easier to keep track of the results, 
## we store them in Tables and then make the plot.

## ---------- PLOT 5.1 RESULTS ---------- ##

## MR - varying the missingness proportion.

## Create the table.
MRTable1m <- matrix(0, 35, 5)
colnames(MRTable1m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable1m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - 90% missing data.
load("MrSimR1_1_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[1:5, 3] <- stderr
MRTable1m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - 80% missing data.
load("MrSimR1_1_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[6:10, 3] <- stderr
MRTable1m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - 65% missing data.
load("MrSimR1_1_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[11:15, 3] <- stderr
MRTable1m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - 50% missing data.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[16:20, 3] <- stderr
MRTable1m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - 35% missing data.
load("MrSimR1_1_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[21:25, 3] <- stderr
MRTable1m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - 80% missing data.
load("MrSimR1_1_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[26:30, 3] <- stderr
MRTable1m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - 10% missing data.
load("MrSimR1_1_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[31:35, 3] <- stderr
MRTable1m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## ---------- PLOT 5.2 RESULTS ---------- ##

## MR - varying the Z-X association.

## Create the table.
MRTable2m <- matrix(0, 30, 5)
colnames(MRTable2m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable2m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.x = 0.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[1:5, 3] <- stderr
MRTable2m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.x = 0.1.
load("MrSimR1_2_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[6:10, 3] <- stderr
MRTable2m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.x = 0.2.
load("MrSimR1_2_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[11:15, 3] <- stderr
MRTable2m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.x = 0.3.
load("MrSimR1_2_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[16:20, 3] <- stderr
MRTable2m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.x = 0.4.
load("MrSimR1_2_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[21:25, 3] <- stderr
MRTable2m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.x = 0.5.
load("MrSimR1_2_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[26:30, 3] <- stderr
MRTable2m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## ---------- PLOT 5.3 RESULTS ---------- ##

## MR - varying the sample size.

## Create the table.
MRTable5m <- matrix(0, 35, 5)
colnames(MRTable5m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable5m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - n = 1000.
load("MrSimR1_5_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[1:5, 3] <- stderr
MRTable5m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - n = 2000.
load("MrSimR1_5_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[6:10, 3] <- stderr
MRTable5m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - n = 5000.
load("MrSimR1_5_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[11:15, 3] <- stderr
MRTable5m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - n = 10000.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[16:20, 3] <- stderr
MRTable5m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - n = 20000.
load("MrSimR1_5_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[21:25, 3] <- stderr
MRTable5m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - n = 50000.
load("MrSimR1_5_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[26:30, 3] <- stderr
MRTable5m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - n = 100000.
load("MrSimR1_5_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[31:35, 3] <- stderr
MRTable5m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## ---------- PLOT 5.4 RESULTS ---------- ##

## MR - varying the Z-Y association.

## Create the table.
MRTable3m <- matrix(0, 30, 5)
colnames(MRTable3m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable3m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.y = 0.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[1:5, 3] <- stderr
MRTable3m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.y = 0.1.
load("MrSimR1_3_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[6:10, 3] <- stderr
MRTable3m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.y = 0.2.
load("MrSimR1_3_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[11:15, 3] <- stderr
MRTable3m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.y = 0.3.
load("MrSimR1_3_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[16:20, 3] <- stderr
MRTable3m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.y = 0.4.
load("MrSimR1_3_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[21:25, 3] <- stderr
MRTable3m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.y = 0.5.
load("MrSimR1_3_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[26:30, 3] <- stderr
MRTable3m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## ---------- GET DATA SUMMARIES ---------- ##

## Extract summaries to be used for plotting
## (i.e. means and confidence interval limits).

## Proportion selected.
Vec.mrmeans1 <- unname(MRTable1m[, 1])
Vec.mrlb1 <- unname(MRTable1m[, 1] - 1.96 * MRTable1m[, 3])
Vec.mrub1 <- unname(MRTable1m[, 1] + 1.96 * MRTable1m[, 3])
Vec.mrmeans1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans1 <- c(NA, Vec.mrmeans1)
Vec.mrlb1 <- c(NA, Vec.mrlb1)
Vec.mrub1 <- c(NA, Vec.mrub1)

## Z-X effect.
Vec.mrmeans2 <- unname(MRTable2m[, 1])
Vec.mrlb2 <- unname(MRTable2m[, 1] - 1.96 * MRTable2m[, 3])
Vec.mrub2 <- unname(MRTable2m[, 1] + 1.96 * MRTable2m[, 3])
Vec.mrmeans2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrlb2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrub2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrmeans2 <- c(NA, Vec.mrmeans2)
Vec.mrlb2 <- c(NA, Vec.mrlb2)
Vec.mrub2 <- c(NA, Vec.mrub2)

## Z-Y effect.
Vec.mrmeans3 <- unname(MRTable3m[, 1])
Vec.mrlb3 <- unname(MRTable3m[, 1] - 1.96 * MRTable3m[, 3])
Vec.mrub3 <- unname(MRTable3m[, 1] + 1.96 * MRTable3m[, 3])
Vec.mrmeans3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans3 <- c(NA, Vec.mrmeans3)
Vec.mrlb3 <- c(NA, Vec.mrlb3)
Vec.mrub3 <- c(NA, Vec.mrub3)

## Sample size.
Vec.mrmeans5 <- unname(MRTable5m[, 1])
Vec.mrlb5 <- unname(MRTable5m[, 1] - 1.96 * MRTable5m[, 3])
Vec.mrub5 <- unname(MRTable5m[, 1] + 1.96 * MRTable5m[, 3])
Vec.mrmeans5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans5 <- c(NA, Vec.mrmeans5)
Vec.mrlb5 <- c(NA, Vec.mrlb5)
Vec.mrub5 <- c(NA, Vec.mrub5)

## ---------- CREATE THE PLOT ---------- ##

## Set up the plot.
pdf(file = "MrEffects.pdf", width = 10, height = 7)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(2, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Proportion selected.
plot(Vec.mrmeans1, xlim = c(1, 37), ylim = c(-0.4, 0.8), pch = 19, col = Vec.cols, xlab = "Proportion Selected", ylab = expression(theta), main = "Proportion Selected", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.cols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-X Effect
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = expression(theta), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.cols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Different sample sizes.
plot(Vec.mrmeans5, xlim = c(1, 37), ylim = c(-0.8, 1.2), pch = 19, col = Vec.cols, xlab = "Sample Size", ylab = expression(theta), main = "Different Sample Sizes", xaxt = "n", cex.lab = 1.2)
#axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1000", "2000", "5000", "10000", "20000", "50000", "100000"))
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1e3", "2e3", "5e3", "1e4", "2e4", "5e4", "1e5"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb5[i], Vec.mrub5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb5[i], Vec.mrlb5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub5[i], Vec.mrub5[i]), col = Vec.cols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-Y Effect
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = expression(theta), main = "Z-Y Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.cols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

##################################################

##########   SUPPLEMENTARY PLOT 6   ##########

## To make it easier to keep track of the results, 
## we store them in Tables and then make the plot.

## ---------- PLOT 6 RESULTS ---------- ##

## MR - varying the value of the causal effect.

## Create the table.
MRTable4m <- matrix(0, 35, 5)
colnames(MRTable4m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable4m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - beta.y = 0.
load("MrSimR1_4_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[1:5, 3] <- stderr
MRTable4m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - beta.y = 0.05.
load("MrSimR1_4_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[6:10, 3] <- stderr
MRTable4m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - beta.y = 0.1.
load("MrSimR1_4_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[11:15, 3] <- stderr
MRTable4m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - beta.y = 0.15.
load("MrSimR1_4_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[16:20, 3] <- stderr
MRTable4m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - beta.y = 0.2.
load("MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[21:25, 3] <- stderr
MRTable4m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - beta.y = 0.25.
load("MrSimR1_4_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[26:30, 3] <- stderr
MRTable4m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - beta.y = 0.3.
load("MrSimR1_4_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(svyipw1.y[, 4]^2 / svyipw1.x[, 3]^2 + svyipw1.y[, 3]^2 * svyipw1.x[, 4]^2 / svyipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(svyipw1.y[, 3] / svyipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                         median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(svyipw1.y[, 3] / svyipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                         sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[31:35, 3] <- stderr
MRTable4m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                         mean( ((svyipw1.y[, 3] / svyipw1.x[, 3]) - 1.96 * (stderr[2])) * ((svyipw1.y[, 3] / svyipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                         mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                         mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                         mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## ---------- GET DATA SUMMARIES ---------- ##

## Extract summaries to be used for plotting
## (i.e. means and confidence interval limits).

## Varying causal effect.
True.thetas <- rep(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), each = 5)
Vec.mrmeans4 <- unname(MRTable4m[, 1]) - True.thetas
Vec.mrlb4 <- unname(MRTable4m[, 1] - 1.96 * MRTable4m[, 3]) - True.thetas
Vec.mrub4 <- unname(MRTable4m[, 1] + 1.96 * MRTable4m[, 3]) - True.thetas
Vec.mrmeans4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans4 <- c(NA, Vec.mrmeans4)
Vec.mrlb4 <- c(NA, Vec.mrlb4)
Vec.mrub4 <- c(NA, Vec.mrub4)
Vec.mrpower4 <- MRTable4m[, 5]
Vec.mrpower4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrpower4 <- c(NA, Vec.mrpower4)

## ---------- CREATE THE PLOT ---------- ##

## Set up the plot.
pdf(file = "MrPowers.pdf", width = 10, height = 4)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(1, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Varying causal effect - Confidence intervals.
plot(Vec.mrmeans4, xlim = c(1, 37), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = expression(theta), ylab = "Bias", main = "Bias", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.cols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## Varying causal effect - Power figures.
plot(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 2], ylim = c(0, 1), xlim = c(0.5, 7.5), pch = 19, type = "b", xlab = expression(theta), ylab = "Power", main = "Empirical Power", xaxt = "n", cex.lab = 1.2, col = colors[1])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 3], pch = 19, type = "b", col = colors[2])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 4], pch = 19, type = "b", col = colors[3])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 5], pch = 19, type = "b", col = colors[4])
axis(side = 1, at = 1:7, labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
legend(x = 0.5, y = 1, legend = c("CCA", "IPW", "Heckman", "TTW"), pch = 19, col = colors[1:4], cex = 0.6)

## Goodbye.
dev.off()

##################################################
