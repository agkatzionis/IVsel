
##########   SIMULATION SUMMARY   ##########

## Here we assess the performance of the new 
## simulations and create tables of results.

## Set working directory.
setwd("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1")

## Load R packages and relevant functions.
load("IVsel_Functions.RData")
library(sampleSelection)
library(ivreg)
library(MendelianRandomization)
library(truncnorm)
library(xtable)

## There are two new observational simulations, 
## four MR-related simulations and an additional 
## short MR simulation with different sample sizes.

## The results of the two observational simulations 
## can be added to the existing figure, while for 
## the MR simulation a new figure can be created 
## (one-sample MR, single instrument for inference).

## But first, summarize results in Tables.

##################################################

##########   OBSERVATIONAL TABLE 1   ##########

## Observational - varying the Z-X association.

## Create the table.
ObsTable1 <- matrix(0, 30, 5)
colnames(ObsTable1) <- c("Mean", "Emp SD", "StdErr", "Cover", "Power")
rownames(ObsTable1) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Store data for beta.x = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
ObsTable1[1:5, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[1:5, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[1:5, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                       mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                       mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                       mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                       mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for beta.x = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_1_results.RData")
ObsTable1[6:10, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[6:10, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[6:10, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[6:10, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[6:10, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for beta.x = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_2_results.RData")
ObsTable1[11:15, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[11:15, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[11:15, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[11:15, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[11:15, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for beta.x = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_3_results.RData")
ObsTable1[16:20, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[16:20, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[16:20, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[16:20, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[16:20, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for beta.x = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_4_results.RData")
ObsTable1[21:25, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[21:25, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[21:25, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[21:25, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[21:25, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for beta.x = 1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_5_results.RData")
ObsTable1[26:30, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable1[26:30, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable1[26:30, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable1[26:30, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable1[26:30, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   OBSERVATIONAL TABLE 2   ##########

## Observational - varying the Z-Y association.

## Create the table.
ObsTable2 <- matrix(0, 30, 5)
colnames(ObsTable2) <- c("Mean", "Emp SD", "StdErr", "Cover", "Power")
rownames(ObsTable2) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Store data for gamma.y = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
ObsTable2[1:5, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[1:5, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[1:5, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[1:5, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                       mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                       mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                       mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                       mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[1:5, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                       mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                       mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                       mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                       mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for gamma.y = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_1_results.RData")
ObsTable2[6:10, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[6:10, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[6:10, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[6:10, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                        mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                        mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                        mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                        mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[6:10, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                        mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                        mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                        mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                        mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for gamma.y = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_2_results.RData")
ObsTable2[11:15, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[11:15, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[11:15, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[11:15, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                         mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                         mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                         mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                         mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[11:15, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                         mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                         mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                         mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                         mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for gamma.y = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_3_results.RData")
ObsTable2[16:20, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[16:20, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[16:20, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[16:20, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                         mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                         mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                         mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                         mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[16:20, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                         mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                         mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                         mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                         mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for gamma.y = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_4_results.RData")
ObsTable2[21:25, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[21:25, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[21:25, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[21:25, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                         mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                         mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                         mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                         mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[21:25, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                         mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                         mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                         mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                         mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

## Store data for gamma.y = 1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_5_results.RData")
ObsTable2[26:30, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
ObsTable2[26:30, 2] <- c(sd(naive[, 3]), sd(ipw1[, 3]), sd(heckman3[, 3]), sd(partial[, 3]), sd(oracle[, 3]))
ObsTable2[26:30, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
ObsTable2[26:30, 4] <- c(mean( (naive[, 3] - 0.1 - 1.96 * naive[, 4]) * (naive[, 3] - 0.1 + 1.96 * naive[, 4]) < 0 ), 
                         mean( (ipw1[, 3] - 0.1 - 1.96 * ipw1[, 4]) * (ipw1[, 3] - 0.1 + 1.96 * ipw1[, 4]) < 0 ), 
                         mean( (heckman3[, 3] - 0.1 - 1.96 * heckman3[, 4]) * (heckman3[, 3] - 0.1 + 1.96 * heckman3[, 4]) < 0 ), 
                         mean( (partial[, 3] - 0.1 - 1.96 * partial[, 4]) * (partial[, 3] - 0.1 + 1.96 * partial[, 4]) < 0 ), 
                         mean( (oracle[, 3] - 0.1 - 1.96 * oracle[, 4]) * (oracle[, 3] - 0.1 + 1.96 * oracle[, 4]) < 0 ) )
ObsTable2[26:30, 5] <- c(mean( (naive[, 3] - 1.96 * naive[, 4]) * (naive[, 3] + 1.96 * naive[, 4]) > 0 ), 
                         mean( (ipw1[, 3] - 1.96 * ipw1[, 4]) * (ipw1[, 3] + 1.96 * ipw1[, 4]) > 0 ), 
                         mean( (heckman3[, 3] - 1.96 * heckman3[, 4]) * (heckman3[, 3] + 1.96 * heckman3[, 4]) > 0 ), 
                         mean( (partial[, 3] - 1.96 * partial[, 4]) * (partial[, 3] + 1.96 * partial[, 4]) > 0 ), 
                         mean( (oracle[, 3] - 1.96 * oracle[, 4]) * (oracle[, 3] + 1.96 * oracle[, 4]) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 1   ##########

## MR - varying the missingness proportion.

## Create the table.
MRTable1 <- matrix(0, 35, 5)
colnames(MRTable1) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable1) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - 90% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[1:5, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[1:5, 3] <- stderr
MRTable1[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - 80% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[6:10, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[6:10, 3] <- stderr
MRTable1[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - 65% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[11:15, 3] <- stderr
MRTable1[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - 50% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[16:20, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[16:20, 3] <- stderr
MRTable1[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - 35% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[21:25, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[21:25, 3] <- stderr
MRTable1[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - 80% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[26:30, 3] <- stderr
MRTable1[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - 10% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1[31:35, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1[31:35, 3] <- stderr
MRTable1[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 2   ##########

## MR - varying the Z-X association.

## Create the table.
MRTable2 <- matrix(0, 30, 5)
colnames(MRTable2) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable2) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.x = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[1:5, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[1:5, 3] <- stderr
MRTable2[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.x = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[6:10, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                       mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[6:10, 3] <- stderr
MRTable2[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.x = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[11:15, 3] <- stderr
MRTable2[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.x = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[16:20, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[16:20, 3] <- stderr
MRTable2[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.x = 0.4.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[21:25, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[21:25, 3] <- stderr
MRTable2[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.x = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2[26:30, 3] <- stderr
MRTable2[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 3   ##########

## MR - varying the Z-Y association.

## Create the table.
MRTable3 <- matrix(0, 30, 5)
colnames(MRTable3) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable3) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.y = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[1:5, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[1:5, 3] <- stderr
MRTable3[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.y = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[6:10, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                       mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[6:10, 3] <- stderr
MRTable3[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.y = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[11:15, 3] <- stderr
MRTable3[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.y = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[16:20, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[16:20, 3] <- stderr
MRTable3[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.y = 0.4.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[21:25, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[21:25, 3] <- stderr
MRTable3[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.y = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3[26:30, 3] <- stderr
MRTable3[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 4   ##########

## MR - varying the value of the causal effect.

## Create the table.
MRTable4 <- matrix(0, 35, 5)
colnames(MRTable4) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable4) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - beta.y = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[1:5, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[1:5, 3] <- stderr
MRTable4[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - beta.y = 0.05.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[6:10, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                       mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[6:10, 3] <- stderr
MRTable4[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - beta.y = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[11:15, 3] <- stderr
MRTable4[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - beta.y = 0.15.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[16:20, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[16:20, 3] <- stderr
MRTable4[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - beta.y = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[21:25, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[21:25, 3] <- stderr
MRTable4[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - beta.y = 0.25.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[26:30, 3] <- stderr
MRTable4[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - beta.y = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4[31:35, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4[31:35, 3] <- stderr
MRTable4[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 5   ##########

## MR - varying the sample size.

## Create the table.
MRTable5 <- matrix(0, 35, 5)
colnames(MRTable5) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable5) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - n = 1000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[1:5, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                      mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[1:5, 3] <- stderr
MRTable5[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - n = 2000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[6:10, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                       mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[6:10, 3] <- stderr
MRTable5[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - n = 5000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[11:15, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[11:15, 3] <- stderr
MRTable5[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - n = 10000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[16:20, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[16:20, 3] <- stderr
MRTable5[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - n = 20000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[21:25, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[21:25, 3] <- stderr
MRTable5[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - n = 50000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[26:30, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[26:30, 3] <- stderr
MRTable5[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - n = 100000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(mean(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            mean(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            mean(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            mean(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            mean(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5[31:35, 1] <- c(mean(naive.y[, 3] / naive.x[, 3]), mean(ipw1.y[, 3] / ipw1.x[, 3]), mean(heckman3.y[, 3] / heckman3.x[, 3]), 
                        mean(partial.y[, 3] / partial.x[, 3]), mean(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5[31:35, 3] <- stderr
MRTable5[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

## Brief assessment of results.
ObsTable1
ObsTable2
MRTable1
MRTable2
MRTable3
MRTable4
MRTable5

## Print all tables in LaTeX format.
xtable(ObsTable1, digits = 3)
xtable(ObsTable2, digits = 3)
xtable(MRTable1, digits = 3)
xtable(MRTable2, digits = 3)
xtable(MRTable3, digits = 3)
xtable(MRTable4, digits = 3)
xtable(MRTable5, digits = 3)

## Worth recreating the MR tables using the median instead of the mode.

##########   MR SIMULATION TABLE 1m   ##########

## MR - varying the missingness proportion.

## Create the table.
MRTable1m <- matrix(0, 35, 5)
colnames(MRTable1m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable1m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - 90% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                      median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[1:5, 3] <- stderr
MRTable1m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - 80% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[6:10, 3] <- stderr
MRTable1m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - 65% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[11:15, 3] <- stderr
MRTable1m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - 50% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[16:20, 3] <- stderr
MRTable1m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - 35% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[21:25, 3] <- stderr
MRTable1m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - 80% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[26:30, 3] <- stderr
MRTable1m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - 10% missing data.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_1_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable1m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable1m[31:35, 3] <- stderr
MRTable1m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable1m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 2m   ##########

## MR - varying the Z-X association.

## Create the table.
MRTable2m <- matrix(0, 30, 5)
colnames(MRTable2m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable2m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.x = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                      median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[1:5, 3] <- stderr
MRTable2m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.x = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[6:10, 3] <- stderr
MRTable2m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.x = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[11:15, 3] <- stderr
MRTable2m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.x = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[16:20, 3] <- stderr
MRTable2m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.x = 0.4.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[21:25, 3] <- stderr
MRTable2m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.x = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_2_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable2m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable2m[26:30, 3] <- stderr
MRTable2m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable2m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 3m   ##########

## MR - varying the Z-Y association.

## Create the table.
MRTable3m <- matrix(0, 30, 5)
colnames(MRTable3m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable3m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 6)

## Scenario 1 - lambda.y = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                      median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[1:5, 3] <- stderr
MRTable3m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - lambda.y = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[6:10, 3] <- stderr
MRTable3m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - lambda.y = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[11:15, 3] <- stderr
MRTable3m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - lambda.y = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[16:20, 3] <- stderr
MRTable3m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - lambda.y = 0.4.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[21:25, 3] <- stderr
MRTable3m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - lambda.y = 0.5.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_3_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable3m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable3m[26:30, 3] <- stderr
MRTable3m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable3m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 4m   ##########

## MR - varying the value of the causal effect.

## Create the table.
MRTable4m <- matrix(0, 35, 5)
colnames(MRTable4m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable4m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - beta.y = 0.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                      median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[1:5, 3] <- stderr
MRTable4m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - beta.y = 0.05.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[6:10, 3] <- stderr
MRTable4m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - beta.y = 0.1.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[11:15, 3] <- stderr
MRTable4m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - beta.y = 0.15.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[16:20, 3] <- stderr
MRTable4m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - beta.y = 0.2.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[21:25, 3] <- stderr
MRTable4m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - beta.y = 0.25.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[26:30, 3] <- stderr
MRTable4m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - beta.y = 0.3.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_4_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable4m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable4m[31:35, 3] <- stderr
MRTable4m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable4m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

##########   MR SIMULATION TABLE 5m   ##########

## MR - varying the sample size.

## Create the table.
MRTable5m <- matrix(0, 35, 5)
colnames(MRTable5m) <- c("Causal", "Emp SD", "StdErr", "Cover", "Power")
rownames(MRTable5m) <- rep(c("CCA", "IPW", "Heckman", "TTW", "Oracle"), times = 7)

## Scenario 1 - n = 1000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_1_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[1:5, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                      median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[1:5, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                      sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[1:5, 3] <- stderr
MRTable5m[1:5, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[1:5, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                      mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                      mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                      mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                      mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 2 - n = 2000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_2_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[6:10, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                       median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[6:10, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                       sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[6:10, 3] <- stderr
MRTable5m[6:10, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[6:10, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                       mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                       mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                       mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                       mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 3 - n = 5000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[11:15, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[11:15, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[11:15, 3] <- stderr
MRTable5m[11:15, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[11:15, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 4 - n = 10000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Mendelian Randomization/MR3_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[16:20, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[16:20, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[16:20, 3] <- stderr
MRTable5m[16:20, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[16:20, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 5 - n = 20000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_4_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[21:25, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[21:25, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[21:25, 3] <- stderr
MRTable5m[21:25, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[21:25, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 6 - n = 50000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_5_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[26:30, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[26:30, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[26:30, 3] <- stderr
MRTable5m[26:30, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[26:30, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

## Scenario 7 - n = 100000.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/MrSimR1_5_XY_6_results.RData")

## Compute second-order standard errors.
stderr <- c(median(sqrt(naive.y[, 4]^2 / naive.x[, 3]^2 + naive.y[, 3]^2 * naive.x[, 4]^2 / naive.x[, 3]^4)), 
            median(sqrt(ipw1.y[, 4]^2 / ipw1.x[, 3]^2 + ipw1.y[, 3]^2 * ipw1.x[, 4]^2 / ipw1.x[, 3]^4)), 
            median(sqrt(heckman3.y[, 4]^2 / heckman3.x[, 3]^2 + heckman3.y[, 3]^2 * heckman3.x[, 4]^2 / heckman3.x[, 3]^4)), 
            median(sqrt(partial.y[, 4]^2 / partial.x[, 3]^2 + partial.y[, 3]^2 * partial.x[, 4]^2 / partial.x[, 3]^4)), 
            median(sqrt(oracle.y[, 4]^2 / oracle.x[, 3]^2 + oracle.y[, 3]^2 * oracle.x[, 4]^2 / oracle.x[, 3]^4)) )

## Store relevant values.
MRTable5m[31:35, 1] <- c(median(naive.y[, 3] / naive.x[, 3]), median(ipw1.y[, 3] / ipw1.x[, 3]), median(heckman3.y[, 3] / heckman3.x[, 3]), 
                        median(partial.y[, 3] / partial.x[, 3]), median(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[31:35, 2] <- c(sd(naive.y[, 3] / naive.x[, 3]), sd(ipw1.y[, 3] / ipw1.x[, 3]), sd(heckman3.y[, 3] / heckman3.x[, 3]), 
                        sd(partial.y[, 3] / partial.x[, 3]), sd(oracle.y[, 3] / oracle.x[, 3]) )
MRTable5m[31:35, 3] <- stderr
MRTable5m[31:35, 4] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 0.2 - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) - 0.2 + 1.96 * (stderr[1])) < 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) - 0.2 + 1.96 * (stderr[2])) < 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) - 0.2 + 1.96 * (stderr[3])) < 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 0.2 - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) - 0.2 + 1.96 * (stderr[4])) < 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) - 0.2 + 1.96 * (stderr[5])) < 0 ) )
MRTable5m[31:35, 5] <- c(mean( ((naive.y[, 3] / naive.x[, 3]) - 1.96 * (stderr[1])) * ((naive.y[, 3] / naive.x[, 3]) + 1.96 * (stderr[1])) > 0 ), 
                        mean( ((ipw1.y[, 3] / ipw1.x[, 3]) - 1.96 * (stderr[2])) * ((ipw1.y[, 3] / ipw1.x[, 3]) + 1.96 * (stderr[2])) > 0 ), 
                        mean( ((heckman3.y[, 3] / heckman3.x[, 3]) - 1.96 * (stderr[3])) * ((heckman3.y[, 3] / heckman3.x[, 3]) + 1.96 * (stderr[3])) > 0 ), 
                        mean( ((partial.y[, 3] / partial.x[, 3]) - 1.96 * (stderr[4])) * ((partial.y[, 3] / partial.x[, 3]) + 1.96 * (stderr[4])) > 0 ), 
                        mean( ((oracle.y[, 3] / oracle.x[, 3]) - 1.96 * (stderr[5])) * ((oracle.y[, 3] / oracle.x[, 3]) + 1.96 * (stderr[5])) > 0 ) )

##################################################

## Brief assessment of results.
cbind(MRTable1[, c(1, 3)], MRTable1m[, c(1, 3)])
cbind(MRTable2[, c(1, 3)], MRTable2m[, c(1, 3)])
cbind(MRTable3[, c(1, 3)], MRTable3m[, c(1, 3)])
cbind(MRTable4[, c(1, 3)], MRTable4m[, c(1, 3)])
cbind(MRTable5[, c(1, 3)], MRTable5m[, c(1, 3)])

## Now plot the results.

##################################################

##########   OBSERVATIONAL PLOT   ##########

## This extends the existing 2x2 observational plot
## by adding varied Z-X and Z-Y effects.

## The colors we will use for plotting are:
colors <- c("red", "brown", "darkgreen", "blue", "black")

## ---------- PLOT OBS1 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means1 <- matrix(0, 5, 6)
rownames(Means1) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means1) <- c("br = 0", "br = 0.1", "br = 0.2", "br = 0.5", "br = 1", "br = 1.5")
Errors1 <- Means1

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim31_0_results.RData")
Means1[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim31_results.RData")
Means1[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim32_results.RData")
Means1[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means1[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim33_results.RData")
Means1[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim34_results.RData")
Means1[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors1[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn them to vectors to make plotting more convenient.
Vec.means1 <- as.vector(rbind(NA, Means1[1:4, ]))
Vec.lb1 <- as.vector(rbind(NA, Means1[1:4, ] - 1.96 * Errors1[1:4, ]))
Vec.ub1 <- as.vector(rbind(NA, Means1[1:4, ] + 1.96 * Errors1[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means1, xlim = c(1, 31), ylim = c(-0.05, 0.25), pch = 19, col = Vec.cols, xlab = "Beta.r", ylab = "Beta.y", main = "X-R Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb1[i], Vec.ub1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb1[i], Vec.lb1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub1[i], Vec.ub1[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT OBS2 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means2 <- matrix(0, 5, 6)
rownames(Means2) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means2) <- c("dr = 0", "dr = 0.1", "dr = 0.2", "dr = 0.5", "dr = 1", "dr = 1.5")
Errors2 <- Means2

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim35_0_results.RData")
Means2[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim35_results.RData")
Means2[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim36_results.RData")
Means2[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means2[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim37_results.RData")
Means2[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim38_results.RData")
Means2[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors2[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means2 and Errors2 to vectors.
Vec.means2 <- as.vector(rbind(NA, Means2[1:4, ]))
Vec.lb2 <- as.vector(rbind(NA, Means2[1:4, ] - 1.96 * Errors2[1:4, ]))
Vec.ub2 <- as.vector(rbind(NA, Means2[1:4, ] + 1.96 * Errors2[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means2, xlim = c(1, 31), ylim = c(-0.05, 0.25), pch = 19, col = Vec.cols, xlab = "Delta.r", ylab = "Beta.y", main = "Y-R Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb2[i], Vec.ub2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb2[i], Vec.lb2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub2[i], Vec.ub2[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT OBS3 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means3 <- matrix(0, 5, 7)
rownames(Means3) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means3) <- c("gr = 0", "gr = 0.08", "gr = 0.2", "gr = 0.27", "gr = 0.4", "gr = 0.6", "gr = 0.95")
Errors3 <- Means3

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim39_results.RData")
Means3[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim40_results.RData")
Means3[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim41_results.RData")
Means3[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim42_results.RData")
Means3[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means3[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim43_results.RData")
Means3[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim44_results.RData")
Means3[, 7] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors3[, 7] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means3 and Errors3 to vectors.
Vec.means3 <- as.vector(rbind(NA, Means3[1:4, ]))
Vec.lb3 <- as.vector(rbind(NA, Means3[1:4, ] - 1.96 * Errors3[1:4, ]))
Vec.ub3 <- as.vector(rbind(NA, Means3[1:4, ] + 1.96 * Errors3[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Do the plotting.
plot(Vec.means3, xlim = c(1, 37), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = "McFadden's R2", ylab = "Beta.y", main = "Instrument Strength", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.1%", "0.5%", "1%", "2%", "5%", "10%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb3[i], Vec.ub3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb3[i], Vec.lb3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub3[i], Vec.ub3[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT OBS4 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means4 <- matrix(0, 5, 7)
rownames(Means4) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means4) <- c("prop = 0.1", "prop = 0.2", "prop = 0.35", "prop = 0.5", "prop = 0.65", "prop = 0.8", "prop = 0.9")
Errors4 <- Means4

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim54_results.RData")
Means4[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim55_results.RData")
Means4[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim56_results.RData")
Means4[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means4[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim57_results.RData")
Means4[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim58_results.RData")
Means4[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim59_results.RData")
Means4[, 7] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors4[, 7] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means4 and Errors4 to vectors.
Vec.means4 <- as.vector(rbind(NA, Means4[1:4, ]))
Vec.lb4 <- as.vector(rbind(NA, Means4[1:4, ] - 1.96 * Errors4[1:4, ]))
Vec.ub4 <- as.vector(rbind(NA, Means4[1:4, ] + 1.96 * Errors4[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Do the plotting.
plot(Vec.means4, xlim = c(1, 37), ylim = c(-0.2, 0.4), pch = 19, col = Vec.cols, xlab = "Prop Selected", ylab = "Beta.y", main = "Proportion Selected", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb4[i], Vec.ub4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb4[i], Vec.lb4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub4[i], Vec.ub4[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT OBS5 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means5 <- matrix(0, 5, 6)
rownames(Means5) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means5) <- c("bx = 0", "bx = 0.1", "bx = 0.2", "bx = 0.3", "bx = 0.5", "bx = 1")
Errors5 <- Means5

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means5[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_1_results.RData")
Means5[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_2_results.RData")
Means5[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_3_results.RData")
Means5[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_4_results.RData")
Means5[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_1_5_results.RData")
Means5[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors5[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means5 and Errors5 to vectors.
Vec.means5 <- as.vector(rbind(NA, Means5[1:4, ]))
Vec.lb5 <- as.vector(rbind(NA, Means5[1:4, ] - 1.96 * Errors5[1:4, ]))
Vec.ub5 <- as.vector(rbind(NA, Means5[1:4, ] + 1.96 * Errors5[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means5, xlim = c(1, 31), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = "Beta.x", ylab = "Beta.y", main = "Z-X Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.5", "1"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb5[i], Vec.ub5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb5[i], Vec.lb5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub5[i], Vec.ub5[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT OBS6 ---------- ##

## Store mean estimates and standard errors from the simulations.
Means6 <- matrix(0, 5, 6)
rownames(Means6) <- c("CC", "IPW", "Heck", "TTW", "Oracle")
colnames(Means6) <- c("gy = 0", "gy = 0.1", "gy = 0.2", "gy = 0.3", "gy = 0.5", "gy = 1")
Errors6 <- Means6

## Load the data and store the appropriate values.
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Observational Case - Final/Sim1_results.RData")
Means6[, 1] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 1] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_1_results.RData")
Means6[, 2] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 2] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_2_results.RData")
Means6[, 3] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 3] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_3_results.RData")
Means6[, 4] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 4] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_4_results.RData")
Means6[, 5] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 5] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))
load("C:/Users/bj20642/OneDrive - University of Bristol/Desktop/Bristol HPC/IVsel - Revision 1/ObsSimR1_2_5_results.RData")
Means6[, 6] <- c(mean(naive[, 3]), mean(ipw1[, 3]), mean(heckman3[, 3]), mean(partial[, 3]), mean(oracle[, 3]))
Errors6[, 6] <- c(mean(naive[, 4]), mean(ipw1[, 4]), mean(heckman3[, 4]), mean(partial[, 4]), mean(oracle[, 4]))

## Turn Means6 and Errors6 to vectors.
Vec.means6 <- as.vector(rbind(NA, Means6[1:4, ]))
Vec.lb6 <- as.vector(rbind(NA, Means6[1:4, ] - 1.96 * Errors6[1:4, ]))
Vec.ub6 <- as.vector(rbind(NA, Means6[1:4, ] + 1.96 * Errors6[1:4, ]))
Vec.cols <- rep(c(NA, colors[1:4]), times = 6)

## Do the plotting.
plot(Vec.means6, xlim = c(1, 31), ylim = c(-0.6, 0.2), pch = 19, col = Vec.cols, xlab = "Gamma.y", ylab = "Beta.y", main = "Z-Y Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.5", "1"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb6[i], Vec.ub6[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb6[i], Vec.lb6[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub6[i], Vec.ub6[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## ---------- OLD PLOT ---------- ##

## We recreate the old plot from the previous files.

## Set up the plot.
pdf(file = "SelectionEffects.pdf", width = 10, height = 7)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(2, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## First, we vary the proportion of individuals with fully observed data.
plot(Vec.means4, xlim = c(1, 37), ylim = c(-0.2, 0.4), pch = 19, col = Vec.cols, xlab = "Proportion Selected", ylab = expression(beta[Y]), main = "Proportion Selected", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
#axis(side = 2, at = c(-2:4 / 10), labels = c("-0.2", "", "0", "", "0.2", "", "0.4"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb4[i], Vec.ub4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb4[i], Vec.lb4[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub4[i], Vec.ub4[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Second, we vary the strength of the covariate-selection association. 
plot(Vec.means1, xlim = c(1, 31), ylim = c(-0.05, 0.25), pch = 19, col = Vec.cols, xlab = expression(beta[R]), ylab = expression(beta[Y]), main = "X-R Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
#axis(side = 2, at = c(-1:5 / 20), labels = c("-0.05", "", "0.05", "", "0.15", "", "0.25"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb1[i], Vec.ub1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb1[i], Vec.lb1[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub1[i], Vec.ub1[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Third, we vary the strength of the instrument-selection association.
## This was already explored in Table 2 earlier and is repeated here.
plot(Vec.means3, xlim = c(1, 37), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = expression(paste(R^2, " Statistic")), ylab = expression(beta[Y]), main = "Instrument Strength", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.1%", "0.5%", "1%", "2%", "5%", "10%"))
#axis(side = 2, at = c(-1:3 / 10), labels = c("-0.1", "0", "0.1", "0.2", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.lb3[i], Vec.ub3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb3[i], Vec.lb3[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub3[i], Vec.ub3[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Fourth, we vary the strength of the outcome-selection association.
plot(Vec.means2, xlim = c(1, 31), ylim = c(-0.05, 0.25), pch = 19, col = Vec.cols, xlab = expression(delta[R]), ylab = expression(beta[Y]), main = "Y-R Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.5", "1", "1.5"))
#axis(side = 2, at = c(-1:5 / 20), labels = c("-0.05", "", "0.05", "", "0.15", "", "0.25"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb2[i], Vec.ub2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb2[i], Vec.lb2[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub2[i], Vec.ub2[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

## ---------- NEW PLOT ---------- ##

## We now plot Z-X and Z-Y effects in a separate figure.

## Set up the plot.
pdf(file = "InstrumentEffects.pdf", width = 10, height = 4)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(1, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## First, plot simulations with a Z-X effect.
plot(Vec.means5, xlim = c(1, 31), ylim = c(-0.1, 0.3), pch = 19, col = Vec.cols, xlab = expression(beta[X]), ylab = expression(beta[Y]), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.5", "1"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.lb5[i], Vec.ub5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.lb5[i], Vec.lb5[i]), col = Vec.cols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.ub5[i], Vec.ub5[i]), col = Vec.cols[i])
}
abline(h = 0.1, col = "gray", lwd = 2, lty = 2)

## Second, plot simulations with a Z-Y effect.
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

##########   MR PLOT   ##########

## We now do a similar plot for the MR analyses.

## We will use the same colors for plotting:
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.mrcols <- rep(c(NA, colors[1:4]), times = 7)

## First, plot each MR run individually.
par(mfrow = c(1, 1))

## ---------- PLOT MR1 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans1 <- unname(MRTable1[, 1])
Vec.mrlb1 <- unname(MRTable1[, 1] - 1.96 * MRTable1[, 3])
Vec.mrub1 <- unname(MRTable1[, 1] + 1.96 * MRTable1[, 3])
Vec.mrmeans1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans1 <- c(NA, Vec.mrmeans1)
Vec.mrlb1 <- c(NA, Vec.mrlb1)
Vec.mrub1 <- c(NA, Vec.mrub1)

## Do the plotting.
plot(Vec.mrmeans1, xlim = c(1, 37), ylim = c(-0.4, 0.8), pch = 19, col = Vec.cols, xlab = "Prob Selected", ylab = "Causal Effect", main = "Proportion Selected", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR2 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans2 <- unname(MRTable2[, 1])
Vec.mrlb2 <- unname(MRTable2[, 1] - 1.96 * MRTable2[, 3])
Vec.mrub2 <- unname(MRTable2[, 1] + 1.96 * MRTable2[, 3])
Vec.mrmeans2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrlb2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrub2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrmeans2 <- c(NA, Vec.mrmeans2)
Vec.mrlb2 <- c(NA, Vec.mrlb2)
Vec.mrub2 <- c(NA, Vec.mrub2)

## Do the plotting.
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = "Causal Effect", main = "Z-X Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR3 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans3 <- unname(MRTable3[, 1])
Vec.mrlb3 <- unname(MRTable3[, 1] - 1.96 * MRTable3[, 3])
Vec.mrub3 <- unname(MRTable3[, 1] + 1.96 * MRTable3[, 3])
Vec.mrmeans3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans3 <- c(NA, Vec.mrmeans3)
Vec.mrlb3 <- c(NA, Vec.mrlb3)
Vec.mrub3 <- c(NA, Vec.mrub3)

## Do the plotting.
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = "Causal Effect", main = "Z-Y Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR4 ---------- ##

## Get means and confidence interval limits.
True.thetas <- rep(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), each = 5)
Vec.mrmeans4 <- unname(MRTable4[, 1]) - True.thetas
Vec.mrlb4 <- unname(MRTable4[, 1] - 1.96 * MRTable4[, 3]) - True.thetas
Vec.mrub4 <- unname(MRTable4[, 1] + 1.96 * MRTable4[, 3]) - True.thetas
Vec.mrmeans4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans4 <- c(NA, Vec.mrmeans4)
Vec.mrlb4 <- c(NA, Vec.mrlb4)
Vec.mrub4 <- c(NA, Vec.mrub4)

## Do the plotting.
plot(Vec.mrmeans4, xlim = c(1, 37), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Theta", ylab = "Bias", main = "Varying Causal Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## ---------- JOINT MR PLOT ---------- ##

## Now put the four plots together.

## Set up the plot.
pdf(file = "MrEffects0.pdf", width = 10, height = 7)

par(mar = c(5.1, 5.1, 4.1, 2.1))
par(mfrow = c(2, 2))
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.cols <- rep(c(NA, colors[1:4]), times = 7)

## Proportion selected.
plot(Vec.mrmeans1, xlim = c(1, 37), ylim = c(-0.4, 0.8), pch = 19, col = Vec.cols, xlab = "Proportion Selected", ylab = expression(theta), main = "Proportion Selected", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-X Effect
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = expression(theta), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Varying causal effect.
plot(Vec.mrmeans4, xlim = c(1, 37), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = expression(theta), ylab = "Bias", main = "Varying Causal Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## Z-Y Effect
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = expression(theta), main = "Z-Y Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

##################################################

##########   MRm PLOT   ##########

## We now do the same using median values, not means.
par(mfrow = c(1, 1))

## We will use the same colors for plotting:
colors <- c("red", "brown", "darkgreen", "blue", "black")
Vec.mrcols <- rep(c(NA, colors[1:4]), times = 7)

## ---------- PLOT MR1 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans1 <- unname(MRTable1m[, 1])
Vec.mrlb1 <- unname(MRTable1m[, 1] - 1.96 * MRTable1m[, 3])
Vec.mrub1 <- unname(MRTable1m[, 1] + 1.96 * MRTable1m[, 3])
Vec.mrmeans1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub1[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans1 <- c(NA, Vec.mrmeans1)
Vec.mrlb1 <- c(NA, Vec.mrlb1)
Vec.mrub1 <- c(NA, Vec.mrub1)

## Do the plotting.
plot(Vec.mrmeans1, xlim = c(1, 37), ylim = c(-0.4, 0.8), pch = 19, col = Vec.cols, xlab = "Prob Selected", ylab = "Causal Effect", main = "Proportion Selected", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("10%", "20%", "35%", "50%", "65%", "80%", "90%"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR2 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans2 <- unname(MRTable2m[, 1])
Vec.mrlb2 <- unname(MRTable2m[, 1] - 1.96 * MRTable2m[, 3])
Vec.mrub2 <- unname(MRTable2m[, 1] + 1.96 * MRTable2m[, 3])
Vec.mrmeans2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrlb2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrub2[c(5, 10, 15, 20, 25, 30)] <- NA
Vec.mrmeans2 <- c(NA, Vec.mrmeans2)
Vec.mrlb2 <- c(NA, Vec.mrlb2)
Vec.mrub2 <- c(NA, Vec.mrub2)

## Do the plotting.
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = "Causal Effect", main = "Z-X Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR3 ---------- ##

## Get means and confidence interval limits.
Vec.mrmeans3 <- unname(MRTable3m[, 1])
Vec.mrlb3 <- unname(MRTable3m[, 1] - 1.96 * MRTable3m[, 3])
Vec.mrub3 <- unname(MRTable3m[, 1] + 1.96 * MRTable3m[, 3])
Vec.mrmeans3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub3[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans3 <- c(NA, Vec.mrmeans3)
Vec.mrlb3 <- c(NA, Vec.mrlb3)
Vec.mrub3 <- c(NA, Vec.mrub3)

## Do the plotting.
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = "Causal Effect", main = "Z-Y Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR4 ---------- ##

## Get means and confidence interval limits.
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

## Do the plotting.
plot(Vec.mrmeans4, xlim = c(1, 37), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Theta", ylab = "Bias", main = "Varying Causal Effect", xaxt = "n")
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## ---------- PLOT MR5 ---------- ##

## We do this for completeness.

## Get means and confidence interval limits.
Vec.mrmeans5 <- unname(MRTable5m[, 1])
Vec.mrlb5 <- unname(MRTable5m[, 1] - 1.96 * MRTable5m[, 3])
Vec.mrub5 <- unname(MRTable5m[, 1] + 1.96 * MRTable5m[, 3])
Vec.mrmeans5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrlb5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrub5[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrmeans5 <- c(NA, Vec.mrmeans5)
Vec.mrlb5 <- c(NA, Vec.mrlb5)
Vec.mrub5 <- c(NA, Vec.mrub5)

## Do the plotting.
plot(Vec.mrmeans5, xlim = c(1, 37), ylim = c(-0.8, 1.2), pch = 19, col = Vec.cols, xlab = "Sample Size", ylab = "Theta", main = "Different Sample Sizes", xaxt = "n")
#axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1000", "2000", "5000", "10000", "20000", "50000", "100000"))
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1e3", "2e3", "5e3", "1e4", "2e4", "5e4", "1e5"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb5[i], Vec.mrub5[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb5[i], Vec.mrlb5[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub5[i], Vec.mrub5[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## ---------- JOINT MR PLOT ---------- ##

## Now put the four plots together.

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
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-X Effect
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = expression(theta), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Varying causal effect.
plot(Vec.mrmeans4, xlim = c(1, 37), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = expression(theta), ylab = "Bias", main = "Varying Causal Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## Z-Y Effect
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = expression(theta), main = "Z-Y Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

## ---------- POWER PLOT ---------- ##

## Just for more information, we provide a scatter-plot with power figures.
Vec.mrpower4 <- MRTable4[, 5]
Vec.mrpower4[c(5, 10, 15, 20, 25, 30, 35)] <- NA
Vec.mrpower4 <- c(NA, Vec.mrpower4)

## Do the plot.
plot(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 2], ylim = c(0, 1), xlim = c(0.5, 7.5), pch = 19, type = "b", xlab = expression(theta), ylab = "Empirical Power", main = "Varying Causal Effect", xaxt = "n", cex.lab = 1.2, col = colors[1])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 3], pch = 19, type = "b", col = colors[2])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 4], pch = 19, type = "b", col = colors[3])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 5], pch = 19, type = "b", col = colors[4])
axis(side = 1, at = 1:7, labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
legend(x = 0.5, y = 1, legend = c("CCA", "IPW", "Heckman", "TTW"), pch = 19, col = colors[1:4])

##################################################

##########   MODIFIED PLOTS   ##########

## To include the simulation with various sample sizes, we modify the plots.

## ---------- NEW JOINT MR PLOT ---------- ##

## Replace the power plot with the sample size one.

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
  lines(c(i, i), c(Vec.mrlb1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb1[i], Vec.mrlb1[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub1[i], Vec.mrub1[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-X Effect
plot(Vec.mrmeans2, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-X Effect", ylab = expression(theta), main = "Z-X Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb2[i], Vec.mrlb2[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub2[i], Vec.mrub2[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Different sample sizes.
plot(Vec.mrmeans5, xlim = c(1, 37), ylim = c(-0.8, 1.2), pch = 19, col = Vec.cols, xlab = "Sample Size", ylab = expression(theta), main = "Different Sample Sizes", xaxt = "n", cex.lab = 1.2)
#axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1000", "2000", "5000", "10000", "20000", "50000", "100000"))
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5, 33.5), labels = c("1e3", "2e3", "5e3", "1e4", "2e4", "5e4", "1e5"))
for (i in 1:36) {
  lines(c(i, i), c(Vec.mrlb5[i], Vec.mrub5[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb5[i], Vec.mrlb5[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub5[i], Vec.mrub5[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Z-Y Effect
plot(Vec.mrmeans3, xlim = c(1, 31), ylim = c(-1, 1), pch = 19, col = Vec.cols, xlab = "Z-Y Effect", ylab = expression(theta), main = "Z-Y Effect", xaxt = "n", cex.lab = 1.2)
axis(side = 1, at = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"))
for (i in 1:30) {
  lines(c(i, i), c(Vec.mrlb3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb3[i], Vec.mrlb3[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub3[i], Vec.mrub3[i]), col = Vec.mrcols[i])
}
abline(h = 0.2, col = "gray", lwd = 2, lty = 2)

## Goodbye.
dev.off()

## ---------- DOUBLE POWER PLOT ---------- ##

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
  lines(c(i, i), c(Vec.mrlb4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrlb4[i], Vec.mrlb4[i]), col = Vec.mrcols[i])
  lines(c(i - 0.2, i + 0.2), c(Vec.mrub4[i], Vec.mrub4[i]), col = Vec.mrcols[i])
}
abline(h = 0, col = "gray", lwd = 2, lty = 2)

## Varying causal effect - Power figures.
plot(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 2], ylim = c(0, 1), xlim = c(0.5, 7.5), pch = 19, type = "b", xlab = expression(theta), ylab = "Power", main = "Empirical Power", xaxt = "n", cex.lab = 1.2, col = colors[1])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 3], pch = 19, type = "b", col = colors[2])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 4], pch = 19, type = "b", col = colors[3])
points(x = 1:7, y = Vec.mrpower4[0:6 * 5 + 5], pch = 19, type = "b", col = colors[4])
axis(side = 1, at = 1:7, labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3"))
legend(x = 0.5, y = 1, legend = c("CCA", "IPW", "Heckman", "TTW"), pch = 19, col = colors[1:4])

## Goodbye.
dev.off()

##################################################
