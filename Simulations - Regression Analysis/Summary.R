
##########   SUMMARY   ##########

## Here we assess the performance of IVsel methods in a set  
## of simulations with regression-type analyses and missing 
## data, and create the Tables included in the paper.

## Set working directory.
setwd( "MY WORKING DIRECTORY" )

## Load R packages and relevant functions.
load("IVsel_Functions.RData")
library(sampleSelection)
library(survey)
library(xtable)

##################################################

## Simulation results are supplied in separate files. 
## Here, we combine those files to create Tables 1 and 2,
## as well as Supplementary Tables 1-4 and Supplementary 
## Figures 1-2 of the paper.


## In each Table we report parameter estimates, the standard
## deviation of estimates, model-based standard errors, 
## empirical coverage/Type I error rates and  (for beta != 0)
## the empirical power to identify a non-null association.

##################################################

##########   TABLE 1   ##########

## This is the main simulations table 
## for regression-type simulations.

## Create the table.
Table1 <- matrix(0, 25, 10)
rownames(Table1) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), 5)
colnames(Table1) <- c("Mean", "Emp SD", "Std Error", "Coverage", "Power", "Mean", "Emp SD", "Std Error", "Coverage", "Power")

## Table 1 contains simulation results from ten different 
## scenarios, which are outlined in the paper.

## Scenario 1: baseline scenario.
load("Reg1_results.RData")
Table1[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 2: null X-Y association
load("Reg2_results.RData")
Table1[1:5, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[1:5, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[1:5, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[1:5, 9] <- c(mean( (naive[, 3] - 0 - 1.96 * naive[, 4]) * (naive[, 3] - 0 + 1.96 * naive[, 4]) < 0 ), 
                mean( (svyipw1[, 3] - 0 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0 + 1.96 * svyipw1[, 4]) < 0 ), 
                mean( (heckman3[, 3] - 0 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0 + 1.96 * heckman3[, 4]) < 0 ), 
                mean( (partial[, 3] - 0 - 1.96 * partial[, 4]) * (partial[, 3] - 0 + 1.96 * partial[, 4]) < 0 ), 
                mean( (oracle[, 3] - 0 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0 + 1.96 * oracle[, 4]) < 0 ) )
Table1[1:5, 10] <- rep(NA, 5)

## Scenario 3: binary instrument.
load("Reg3_results.RData")
Table1[6:10, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[6:10, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[6:10, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[6:10, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[6:10, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                     mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                     mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                     mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                     mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 4: binary covariate.
load("Reg4_results.RData")
Table1[6:10, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[6:10, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[6:10, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[6:10, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[6:10, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 5: binary outcome (logistic regression).
load("Reg5_results.RData")
Table1[11:15, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman4[, 3]), mean(full[, 3]), mean(oracle[, 3]))
Table1[11:15, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman4[, 3]), sd(full[, 3]), sd(oracle[, 3]))
Table1[11:15, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman4[, 4]), mean(full[, 4], na.rm = TRUE), mean(oracle[, 4]))
Table1[11:15, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman4[, 3] - 0.1 - 1.96 * heckman4[, 4]) * (heckman4[, 3] - 0.1 + 1.96 * heckman4[, 4]) < 0 ), 
                     mean( (full[, 3] - 0.1 - 1.96 * full[, 4]) * (full[, 3] - 0.1 + 1.96 * full[, 4]) < 0, na.rm = TRUE ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[11:15, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                     mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                     mean( (heckman4[, 3] - 1.96 * heckman4[, 4]) * (heckman4[, 3] + 1.96 * heckman4[, 4]) > 0 ), 
                     mean( (full[, 3] - 1.96 * full[, 4]) * (full[, 3] + 1.96 * full[, 4]) > 0, na.rm = TRUE ), 
                     mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 6: discrete outcome (Poisson regression).
load("Reg6_results.RData")
Table1[11:15, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[11:15, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[11:15, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[11:15, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[11:15, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 7: no X-R effect.
load("Reg7_results.RData")
Table1[16:20, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[16:20, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[16:20, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[16:20, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[16:20, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                     mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                     mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                     mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                     mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 8: no Y-R effect (MAR data).
load("Reg8_results.RData")
Table1[16:20, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[16:20, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[16:20, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[16:20, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[16:20, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 9: instrument affects covariate.
load("Reg9_results.RData")
Table1[21:25, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[21:25, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[21:25, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[21:25, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[21:25, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                     mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                     mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                     mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                     mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 10: instrument affects outcome (IV violation).
load("Reg10_results.RData")
Table1[21:25, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table1[21:25, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table1[21:25, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table1[21:25, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table1[21:25, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   TABLE 2   ##########

## Table 2 explores the performance of IVsel methods
## with varying degrees of instrument strength.

## We quantify instrument strength according to the R^2 
## statistic for the variance in selection explained by 
## the instrument. Since selection is a binary variable,
## we use McFadden's pseudo-R^2 statistic (which is only 
## an approximation of the true variation explained).

## Create the table.
Table2 <- matrix(0, 15, 10)
rownames(Table2) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), 3)
colnames(Table2) <- c("Mean", "Emp SD", "Std Error", "Coverage", "Power", "Mean", "Emp SD", "Std Error", "Coverage", "Power")

## Six scenarios, corresponding to R^2 values of 
## 0.1%, 0.5%, 1%, 2%, 5%, 10%. The 2% simulation 
## is the baseline scenario of Table 1, so we simply
## copy-paste results from that scenario.

## Proportion of variation explained: R^2 = 0.1%.
load("Reg11_results.RData")
Table2[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                    mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                    mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                    mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                    mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                    mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                    mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                    mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                    mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Proportion of variation explained: R^2 = 0.1%.
load("Reg12_results.RData")
Table2[1:5, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[1:5, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[1:5, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[1:5, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                    mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                    mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                    mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                    mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[1:5, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                    mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                    mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                    mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                    mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )


## Proportion of variation explained: R^2 = 0.1%.
load("Reg13_results.RData")
Table2[6:10, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[6:10, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[6:10, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[6:10, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[6:10, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                     mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                     mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                     mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                     mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Proportion of variation explained: R^2 = 0.1%.
load("Reg1_results.RData")
Table2[6:10, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[6:10, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[6:10, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[6:10, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                     mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                     mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                     mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                     mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[6:10, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Proportion of variation explained: R^2 = 0.1%.
load("Reg14_results.RData")
Table2[11:15, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[11:15, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[11:15, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[11:15, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                      mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                      mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                      mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                      mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[11:15, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Proportion of variation explained: R^2 = 0.1%.
load("Reg15_results.RData")
Table2[11:15, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Table2[11:15, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
Table2[11:15, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
Table2[11:15, 9] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                      mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                      mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                      mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                      mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
Table2[11:15, 10] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   SUPPLEMENTARY TABLE 1   ##########

## Supplementary Table 1 reproduces the analyses
## of Table 1 but for a null X-Y association.

## Create the table.
SuppTable1 <- matrix(0, 20, 8)
rownames(SuppTable1) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), 4)
colnames(SuppTable1) <- c("Mean", "Emp SD", "Std Error", "Type I", "Mean", "Emp SD", "Std Error", "Type 1")

## There are eight scenarios (we exclude the
## baseline scenario since null results for 
## that scenario were reported in Table 1).

## Scenario 1: binary instrument.
load("RegS1_results.RData")
SuppTable1[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[1:5, 4] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                    mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                    mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                    mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                    mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 2: binary covariate.
load("RegS2_results.RData")
SuppTable1[1:5, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[1:5, 6] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[1:5, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[1:5, 8] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 3: binary outcome (logistic regression).
load("RegS3_results.RData")
SuppTable1[6:10, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman4[, 3]), mean(full[, 3]), mean(oracle[, 3]))
SuppTable1[6:10, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman4[, 3]), sd(full[, 3]), sd(oracle[, 3]))
SuppTable1[6:10, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman4[, 4]), mean(full[, 4], na.rm = TRUE), mean(oracle[, 4]))
SuppTable1[6:10, 4] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                          mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                          mean( (heckman4[, 3] - 1.96 * heckman4[, 4]) * (heckman4[, 3] + 1.96 * heckman4[, 4]) > 0 ), 
                          mean( (full[, 3] - 1.96 * full[, 4]) * (full[, 3] + 1.96 * full[, 4]) > 0, na.rm = TRUE ), 
                          mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 4: discrete outcome (Poisson regression).
load("RegS4_results.RData")
SuppTable1[6:10, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[6:10, 6] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[6:10, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[6:10, 8] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 5: no X-R effect.
load("RegS5_results.RData")
SuppTable1[11:15, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[11:15, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[11:15, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[11:15, 4] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                          mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                          mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                          mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                          mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 6: no Y-R effect (MAR data).
load("RegS6_results.RData")
SuppTable1[11:15, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[11:15, 6] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[11:15, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[11:15, 8] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 7: instrument affects covariate.
load("RegS7_results.RData")
SuppTable1[16:20, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[16:20, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[16:20, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[16:20, 4] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                      mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                      mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                      mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                      mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 8: instrument affects outcome (IV violation).
load("RegS8_results.RData")
SuppTable1[16:20, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable1[16:20, 6] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable1[16:20, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable1[16:20, 8] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   SUPPLEMENTARY PLOT 1   ##########

## This plots the performance of IVsel methods for
## varying values of the simulation parameters. We 
## do a single plot instead of four separate tables
## for ease of presentation.

## Set up the plot.
pdf(file = "SelectionEffects.pdf", width = 10, height = 7)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(2, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")


## First, we vary the proportion of individuals with fully observed data.

## Store estimates and standard errors from simulations.
Means1 <- matrix(0, 5, 7)
rownames(Means1) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means1) <- c("prop = 0.1", "prop = 0.2", "prop = 0.35", "prop = 0.5", "prop = 0.65", "prop = 0.8", "prop = 0.9")
Errors1 <- Means1

## Load the data and store the appropriate values.
load("RegS9_results.RData")
Means1[, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 1] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS10_results.RData")
Means1[, 2] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 2] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS11_results.RData")
Means1[, 3] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg1_results.RData")
Means1[, 4] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 4] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS12_results.RData")
Means1[, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 5] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS13_results.RData")
Means1[, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 6] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS14_results.RData")
Means1[, 7] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means1 and Errors1 to vectors.
Vec.means1 <- as.vector(rbind(NA, Means1[1:4, ]))
Vec.lb1 <- as.vector(rbind(NA, Means1[1:4, ] - 1.96 * Errors1[1:4, ]))
Vec.ub1 <- as.vector(rbind(NA, Means1[1:4, ] + 1.96 * Errors1[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Do the plotting.
plot(Vec.means1, xlim = c(1, 37), ylim = c(-0.2, 0.4), pch = 19, col = Vec.cols, xlab = "Proportion Selected", ylab = expression(beta[y]), main = "Proportion Selected", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb1[i], Vec.ub1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb1[i], Vec.lb1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub1[i], Vec.ub1[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)


## Second, we vary the strength of the covariate-selection association. 

## Store estimates and standard errors from simulations.
Means2 <- matrix(0, 5, 6)
rownames(Means2) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means2) <- c("br = 0", "br = 0.1", "br = 0.2", "br = 0.5", "br = 1", "br = 1.5")
Errors2 <- Means2

## Load the data and store the appropriate values.
load("RegS15_results.RData")
Means2[, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 1] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS16_results.RData")
Means2[, 2] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 2] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS17_results.RData")
Means2[, 3] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg1_results.RData")
Means2[, 4] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 4] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS18_results.RData")
Means2[, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 5] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS19_results.RData")
Means2[, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 6] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means2 and Errors2 to vectors.
Vec.means2 <- as.vector(rbind(NA, Means2[1:4, ]))
Vec.lb2 <- as.vector(rbind(NA, Means2[1:4, ] - 1.96 * Errors2[1:4, ]))
Vec.ub2 <- as.vector(rbind(NA, Means2[1:4, ] + 1.96 * Errors2[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means2, xlim = c(1, 31), ylim = c(-0.06, 0.25), pch = 19, col = Vec.cols, xlab = expression(beta[r]), ylab = expression(beta[y]), main = "X-R Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb2[i], Vec.ub2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb2[i], Vec.lb2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub2[i], Vec.ub2[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)


## Third, we vary the strength of the instrument-selection association.
## This was already explored in Table 2 earlier and is repeated here.

## Store estimates and standard errors from simulations.
Means3 <- matrix(0, 5, 7)
rownames(Means3) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means3) <- c("gr = 0", "gr = 0.08", "gr = 0.2", "gr = 0.27", "gr = 0.4", "gr = 0.6", "gr = 0.95")
Errors3 <- Means3

## Load the data and store the appropriate values.
load("RegS20_results.RData")
Means3[, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 1] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg11_results.RData")
Means3[, 2] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 2] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg12_results.RData")
Means3[, 3] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg13_results.RData")
Means3[, 4] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 4] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg1_results.RData")
Means3[, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 5] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg14_results.RData")
Means3[, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 6] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg15_results.RData")
Means3[, 7] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 7] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means3 and Errors3 to vectors.
Vec.means3 <- as.vector(rbind(NA, Means3[1:4, ]))
Vec.lb3 <- as.vector(rbind(NA, Means3[1:4, ] - 1.96 * Errors3[1:4, ]))
Vec.ub3 <- as.vector(rbind(NA, Means3[1:4, ] + 1.96 * Errors3[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Do the plotting.
plot(Vec.means3, xlim = c(1, 37), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = expression(paste(R^2, " Statistic")), ylab = expression(beta[y]), main = "Instrument Strength", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.1%", "0.5%", "1%", "2%", "5%", "10%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb3[i], Vec.ub3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb3[i], Vec.lb3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub3[i], Vec.ub3[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)


## Fourth, we vary the strength of the outcome-selection association.

## Store estimates and standard errors from simulations.
Means4 <- matrix(0, 5, 6)
rownames(Means4) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means4) <- c("dr = 0", "dr = 0.1", "dr = 0.2", "dr = 0.5", "dr = 1", "dr = 1.5")
Errors4 <- Means4

## Load the data and store the appropriate values.
load("RegS21_results.RData")
Means4[, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 1] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS22_results.RData")
Means4[, 2] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 2] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS23_results.RData")
Means4[, 3] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("Reg1_results.RData")
Means4[, 4] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 4] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS24_results.RData")
Means4[, 5] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 5] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("RegS25_results.RData")
Means4[, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 6] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means4 and Errors4 to vectors.
Vec.means4 <- as.vector(rbind(NA, Means4[1:4, ]))
Vec.lb4 <- as.vector(rbind(NA, Means4[1:4, ] - 1.96 * Errors4[1:4, ]))
Vec.ub4 <- as.vector(rbind(NA, Means4[1:4, ] + 1.96 * Errors4[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means4, xlim = c(1, 31), ylim = c(-0.05, 0.25), pch = 19, col = Vec.cols, xlab = expression(delta[r]), ylab = expression(beta[y]), main = "Y-R Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb4[i], Vec.ub4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb4[i], Vec.lb4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub4[i], Vec.ub4[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

##################################################

##########   SUPPLEMENTARY PLOT 2   ##########

## This plots the performance of IVsel methods for
## varying values of the Z-X and Z-Y effects.

## Set up the plot.
pdf(file = "SelectionEffects.pdf", width = 10, height = 7)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(1, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")


## First, we vary the Z-X effect.

## Store mean estimates and standard errors from the simulations.
Means5 <- matrix(0, 5, 6)
rownames(Means5) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means5) <- c("bx = 0", "bx = 0.1", "bx = 0.2", "bx = 0.3", "bx = 0.5", "bx = 1")
Errors5 <- Means5

## Load the data and store the appropriate values.
load("Reg1_results.RData")
Means5[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_1_1_results.RData")
Means5[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_1_2_results.RData")
Means5[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_1_3_results.RData")
Means5[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_1_4_results.RData")
Means5[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_1_5_results.RData")
Means5[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means5 and Errors5 to vectors.
Vec.means5 <- as.vector(rbind(NA, Means5[1:4, ]))
Vec.lb5 <- as.vector(rbind(NA, Means5[1:4, ] - 1.96 * Errors5[1:4, ]))
Vec.ub5 <- as.vector(rbind(NA, Means5[1:4, ] + 1.96 * Errors5[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means5, xlim = c(1, 31), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = expression(beta[X]), ylab = expression(beta[Y]), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.5", "1"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb5[i], Vec.ub5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb5[i], Vec.lb5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub5[i], Vec.ub5[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)


## Second, we vary the Z-Y effect.

## Store mean estimates and standard errors from the simulations.
Means6 <- matrix(0, 5, 6)
rownames(Means6) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means6) <- c("gy = 0", "gy = 0.1", "gy = 0.2", "gy = 0.3", "gy = 0.5", "gy = 1")
Errors6 <- Means6

## Load the data and store the appropriate values.
load("Reg1_results.RData")
Means6[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_2_1_results.RData")
Means6[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_2_2_results.RData")
Means6[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_2_3_results.RData")
Means6[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_2_4_results.RData")
Means6[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("ObsSimR1_2_5_results.RData")
Means6[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means6 and Errors6 to vectors.
Vec.means6 <- as.vector(rbind(NA, Means6[1:4, ]))
Vec.lb6 <- as.vector(rbind(NA, Means6[1:4, ] - 1.96 * Errors6[1:4, ]))
Vec.ub6 <- as.vector(rbind(NA, Means6[1:4, ] + 1.96 * Errors6[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means6, xlim = c(1, 31), ylim = c(-0.6, 0.2), pch = 19, col = Vec.cols, xlab = expression(gamma[Y]), ylab = expression(beta[Y]), main = "Z-Y Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.5", "1"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb6[i], Vec.ub6[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb6[i], Vec.lb6[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub6[i], Vec.ub6[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

##################################################

##########   SUPPLEMENTARY TABLE 2   ##########

## Supplementary Table 2 explores the performance of 
## IVsel methods with a weak binary instrument for selection.

## Create the table.
SuppTable2 <- matrix(0, 5, 9)
rownames(SuppTable2) <- c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle")
colnames(SuppTable2) <- c("Mean", "Emp SD", "Std Error", "Coverage", "Power", "Mean", "Emp SD", "Std Error", "Type 1")

## Two scenarios: positive or null X-Y association.

## Scenario 1: positive X-Y association.
load("RegS26_results.RData")
SuppTable2[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable2[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable2[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable2[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                    mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                    mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                    mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                    mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
SuppTable2[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                    mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                    mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                    mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                    mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 2: null X-Y association.
load("RegS27_results.RData")
SuppTable2[1:5, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable2[1:5, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable2[1:5, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable2[1:5, 9] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   SUPPLEMENTARY TABLE 3   ##########

## Supplementary Table 3 explores the performance of 
## IVsel methods when the outcome is confounded with 
## selection instead of directly causing it.

## Create the table.
SuppTable3 <- matrix(0, 5, 9)
rownames(SuppTable3) <- c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle")
colnames(SuppTable3) <- c("Mean", "Emp SD", "Std Error", "Coverage", "Power", "Mean", "Emp SD", "Std Error", "Type 1")

## Two scenarios: positive or null X-Y association.

## Scenario 1: positive X-Y association.
load("RegS28_results.RData")
SuppTable3[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable3[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable3[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable3[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
SuppTable3[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 2: null X-Y association.
load("RegS29_results.RData")
SuppTable3[1:5, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable3[1:5, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable3[1:5, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable3[1:5, 9] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   SUPPLEMENTARY TABLE 4   ##########

## Finally, Supplementary Table 4 explores violations 
## of the normality assumption about the outcome.

## Create the table.
SuppTable4 <- matrix(0, 15, 9)
rownames(SuppTable4) <- rep(c("Complete-case", "IPW", "Heckman", "TTW-Partial", "Oracle"), 3)
colnames(SuppTable4) <- c("Mean", "Emp SD", "Std Error", "Coverage", "Power", "Mean", "Emp SD", "Std Error", "Type 1")

## Six scenarios: t(4), log-normal or bimodal outcome 
## distribution, all implemented for either a positive
## or a null X-Y association.

## Scenario 1: t(4) outcome, positive X-Y association.
load("RegS30_results.RData")
SuppTable4[1:5, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[1:5, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[1:5, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
SuppTable4[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 2: t(4) outcome, null X-Y association.
load("RegS31_results.RData")
SuppTable4[1:5, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[1:5, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[1:5, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[1:5, 9] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 3: log-normal outcome, positive X-Y association.
load("RegS32_results.RData")
SuppTable4[6:10, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[6:10, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[6:10, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[6:10, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
SuppTable4[6:10, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 4: log-normal outcome, null X-Y association.
load("RegS33_results.RData")
SuppTable4[6:10, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[6:10, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[6:10, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[6:10, 9] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 5: bimodal outcome, positive X-Y association.
load("RegS34_results.RData")
SuppTable4[11:15, 1] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[11:15, 2] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[11:15, 3] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[11:15, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (svyipw1[, 3] - 0.1 - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] - 0.1 + 1.96 * svyipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
SuppTable4[11:15, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Scenario 6: bimodal outcome, null X-Y association.
load("RegS35_results.RData")
SuppTable4[11:15, 6] <- c(mean(naive[, 3]), mean(svyipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
SuppTable4[11:15, 7] <- c(sd(naive[, 3]), sd(svyipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
SuppTable4[11:15, 8] <- c(mean(naive[, 4]), mean(svyipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
SuppTable4[11:15, 9] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (svyipw1[, 3] - 1.96 * svyipw1[, 4]) * (svyipw1[, 3] + 1.96 * svyipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

## Print all Tables.
xtable(Table1, digits = 3)
xtable(Table2, digits = 3)
xtable(SuppTable1, digits = 3)
xtable(SuppTable2, digits = 3)
xtable(SuppTable3, digits = 3)
xtable(SuppTable4, digits = 3)

##################################################
