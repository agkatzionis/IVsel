
##########   ALSPAC MR ANALYSIS   ##########

## This script implements a Mendelian randomization 
## analysis of the ALSPAC dataset given to us under
## application B3838. The aim is to assess if BMI 
## is a cause of smoking, self-harm and depression, as
## well as to assess whether selection bias pertinent
## in such analyses can be removed using IVsel methods.
## MR analyses are adjusted for sex and maternal traits. 

## For IVsel methods, we use four instruments for
## selection: one is the randomized trial for ALSPAC 
## participation conducted by Bray et al. (2017), 
## another consisting of variants associated with overall 
## children participation in Taylor et al. (2018) and
## two more PRS derived from optional participation in 
## UKB data (FFQ and MHQ, Tyrell et al. 2021).

##########   ALSPAC SETUP   ##########

## Set working directory, load commands and R packages.
options(digits = 4)
setwd( "[MY WORKING DIRECTORY]" )
load("IVsel_Functions.RData")
library(sampleSelection)
library(estimatr)
library(xtable)


## Load the ALSPAC data.
alspac <- read.csv("Gkatzionis_12Nov2021.csv")

## Load the file that links phenotypic data IDs with genetic data IDs.
linker <- read.csv("B3838_datasetids.csv")

## Load the PRS for BMI, created using variants in the Pulit et al. (2017) GWAS.
bmiscore <- read.table("bmi_scores2.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

## Load genetic PRS to be used as instruments for selection.
prs_h4 <- read.table("Part_prs4_score.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
prs_ffq <- read.table("ffq_score.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
prs_mhq <- read.table("mhq_score.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

## Size of each dataset.
dim(alspac)
dim(linker)
dim(bmiscore)
dim(prs_h4)
dim(prs_ffq)
dim(prs_mhq)


## We need to link the variable cidB3838 between alspac and linker 
## and the variable IID in the PRS files with gi_hrc_g0m_g1 in linker.
linker <- linker[order(linker[, 3]), ]

## Common individuals.
cids <- which(linker$cidB3838 %in% alspac$cidB3838)
cids2 <- which(alspac$cidB3838 %in% linker$cidB3838)
length(cids)
length(cids2)

## ALSPAC contains double entries - find them.
which(diff(alspac$cidB3838) == 0)
cbind(alspac[which(alspac$qlet == "B") - 1, 1:2], alspac[which(alspac$qlet == "B"), 1:2])
## They represent twins, labelled A, B respectively.

## Add linker IDs to ALSPAC for individuals present in both.
alspac <- cbind(alspac, 0, 0)   ## The zeroes are at columns 62 and 63.
colnames(alspac)[c(62, 63)] <- colnames(linker)[c(1, 2)]
for (i in 1:nrow(alspac)) alspac[i, c(62, 63)] <- linker[which(linker$cidB3838 == alspac$cidB3838[i]), c(1, 2)]

## Now link ALSPAC with the Polygenic scores.
## First, create a similar variable in ALSPAC.
alspac <- cbind(alspac, 0)   ## This is column 64.
colnames(alspac)[64] <- "HRC_plus_qlet"
for (i in 1:nrow(alspac)) alspac[i, 64] <- paste(alspac[i, 63], alspac[i, 2], sep = "")

## Test whether this can be matched with polygenic scores.
c(sum(alspac$HRC_plus_qlet %in% bmiscore[, 1]), sum(bmiscore[, 1] %in% alspac$HRC_plus_qlet))
c(sum(alspac$HRC_plus_qlet %in% prs_h4[, 1]), sum(prs_h4[, 1] %in% alspac$HRC_plus_qlet))
c(sum(alspac$HRC_plus_qlet %in% prs_ffq[, 1]), sum(prs_ffq[, 1] %in% alspac$HRC_plus_qlet))
c(sum(alspac$HRC_plus_qlet %in% prs_mhq[, 1]), sum(prs_mhq[, 1] %in% alspac$HRC_plus_qlet))

## Now add the BMI risk score and the three PRS to the ALSPAC dataset.
alspac <- cbind(alspac, NA, NA, NA, NA)   ## These are columns 65-68.
colnames(alspac)[65:68] <- c("BMI_grs", "PRS1", "PRS2", "PRS3")
for (i in 1:nrow(alspac)) {
  if (alspac$HRC_plus_qlet[i] %in% bmiscore[, 1]) alspac[i, 65] <- bmiscore[which(bmiscore[, 1] == alspac$HRC_plus_qlet[i]), 6]
  if (alspac$HRC_plus_qlet[i] %in% prs_h4[, 1]) alspac[i, 66] <- prs_h4[which(prs_h4[, 1] == alspac$HRC_plus_qlet[i]), 6]
  if (alspac$HRC_plus_qlet[i] %in% prs_ffq[, 1]) alspac[i, 67] <- prs_ffq[which(prs_ffq[, 1] == alspac$HRC_plus_qlet[i]), 6]
  if (alspac$HRC_plus_qlet[i] %in% prs_mhq[, 1]) alspac[i, 68] <- prs_mhq[which(prs_mhq[, 1] == alspac$HRC_plus_qlet[i]), 6]
}

## Now leave out individuals with no genetic data.
alspac2 <- alspac[which(!(is.na(alspac$BMI_grs))), ]
dim(alspac2)

##################################################

##########   EXTRACT VARIABLES   ##########

## Store the analysis variables as separate vectors.
## We mostly use the CCU data as they seem to come from  
## the same source as the RCT results.
bmi <- alspac2$FJMR022a
rct <- alspac2$CCUarm
participant <- alspac2$CCU0006
smoking.ever <- alspac2$CCU3000
smoking.freq <- alspac2$CCU3013
harm <- alspac2$CCU2040
depression.all <- alspac2$FJCI1003
birthweight <- alspac2$kz030
sex <- alspac2$kz021
mother.age <- alspac2$mz028b
mother.smk <- alspac2$b663
mother.educ <- alspac2$c645a

## Smoking, Number of cigarettes and Self-Harm are from
## the "It's all about you" stage (20+ years).
## BMI measurements come from TF4 (18 years).

## Some of these have missing values, remove them first.
bmi[bmi < 0] <- NA
rct[rct < 0] <- NA; rct[rct == 2] <- 0   ## Not clear. I guess 1: email + post, 0: email.
participant[participant < 0] <- NA; participant[participant == 2] <- 0   ## 1: yes, 0: no.
smoking.ever[smoking.ever < 0] <- NA; smoking.ever[smoking.ever == 2] <- 0   ## 1: yes, 0: no.
smoking.freq[smoking.freq == -2] <- 0; smoking.freq[smoking.freq == -3] <- 0; smoking.freq[smoking.freq < 0] <- NA
harm[harm < 0] <- NA; harm[harm == 2] <- 0   ## 1: yes, 0: no.
depression.all[depression.all < 0] <- NA   ## 1: yes, 0: no.
birthweight[birthweight < 0] <- NA
sex[sex < 0] <- NA; sex[sex == 2] <- 0   ## 1: male, 0: female
mother.age[mother.age < 0] <- NA
mother.smk[mother.smk < 0] <- NA; mother.smk[mother.smk == 1] <- 0; mother.smk[mother.smk > 1] <- 1   ##   ## 1: yes, 0: no.
mother.educ[mother.educ < 0] <- NA   ## Categorical, 6 categories.

## Identify individuals with missing covariate values.
sum(!(is.na(sex)) & !(is.na(mother.age)) & !(is.na(mother.smk)) & !(is.na(mother.educ)) & !(is.na(rct)))   ## 3175
sum(!(is.na(sex)) & !(is.na(mother.age)) & !(is.na(mother.smk)) & !(is.na(mother.educ)))   ## 4971
ind1 <- which(!(is.na(sex)) & !(is.na(mother.age)) & !(is.na(mother.smk)) & !(is.na(mother.educ)) & !(is.na(rct)) == TRUE)
ind2 <- which(!(is.na(sex)) & !(is.na(mother.age)) & !(is.na(mother.smk)) & !(is.na(mother.educ)) == TRUE)

## Group together the adjustment variables.
covariates <- cbind(sex, mother.age, mother.smk, mother.educ)

## Extract the BMI genetic risk score.
grs <- alspac2$BMI_grs

## Extract the three PRS scores for participation.
prs1 <- alspac2$PRS1
prs2 <- alspac2$PRS2
prs3 <- alspac2$PRS3

## Standardize the risk scores for stability.
grs.std <- (grs - mean(grs, na.rm = TRUE)) / sd(grs, na.rm = TRUE)
prs1.std <- (prs1 - mean(prs1, na.rm = TRUE)) / sd(prs1, na.rm = TRUE)
prs2.std <- (prs2 - mean(prs2, na.rm = TRUE)) / sd(prs2, na.rm = TRUE)
prs3.std <- (prs3 - mean(prs3, na.rm = TRUE)) / sd(prs3, na.rm = TRUE)

##################################################

## Now implement the MR analysis, using either the RCT 
## or the genetic risk scores as instruments for 
## participation. We also fit CCA and IPW for comparison.

## This is done for each trait separately.

##################################################

##########   BODY MASS INDEX   ##########

## For the MR analysis, we need to estimate the GRS-BMI 
## association. We adjust for sex and the three maternal traits.
## BMI is normally distributed so we use linear regression.

## Subset BMI as necessary.
bmi1 <- bmi[ind1]
bmi2 <- bmi[ind2]
R.1 <- !(is.na(bmi1))
R.2 <- !(is.na(bmi2))

## Proportion of missing data.
sum(is.na(bmi1)) / length(bmi1)
sum(is.na(bmi2)) / length(bmi2)

## Complete-case analysis.
bmi.cc <- summary(lm(bmi2 ~ grs.std[ind2] + covariates[ind2, ]))

## Inverse Probability Weighting.
bmi.ipw <- ipw.linear(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = bmi2, R = R.2)$est

## Heckman's selection model.
bmi.heck1 <- summary(selection(R.1 ~ grs.std[ind1] + covariates[ind1, ] + rct[ind1], bmi1 ~ grs.std[ind1] + covariates[ind1, ], method = "2step"))
bmi.heck2 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs1.std[ind2], bmi2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))
bmi.heck3 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs2.std[ind2], bmi2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))
bmi.heck4 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs3.std[ind2], bmi2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))

## Tchetgen Tchetgen and Wirth.
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1])
bmi.ttw.m1 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = bmi1, alpha.hat = partial.opt21$par, hessian = TRUE)
bmi.ttw.s1 <- sqrt(diag(solve(- bmi.ttw.m1$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2])
bmi.ttw.m2 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = bmi2, alpha.hat = partial.opt21$par, hessian = TRUE)
bmi.ttw.s2 <- sqrt(diag(solve(- bmi.ttw.m2$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2])
bmi.ttw.m3 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = bmi2, alpha.hat = partial.opt21$par, hessian = TRUE)
bmi.ttw.s3 <- sqrt(diag(solve(- bmi.ttw.m3$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2])
bmi.ttw.m4 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = bmi2, alpha.hat = partial.opt21$par, hessian = TRUE)
bmi.ttw.s4 <- sqrt(diag(solve(- bmi.ttw.m4$hessian)))

## Assess standardized estimates.
bmi.cc$coefficients
bmi.ipw
bmi.heck1
bmi.heck2
bmi.heck3
bmi.heck4
cbind(bmi.ttw.m1$par, bmi.ttw.s1)
cbind(bmi.ttw.m2$par, bmi.ttw.s2)
cbind(bmi.ttw.m3$par, bmi.ttw.s3)
cbind(bmi.ttw.m4$par, bmi.ttw.s4)

##################################################

##########   SMOKING STATUS   ##########

## Smoking status is binary and will be treated as such.

## Subset BMI as necessary.
smk1 <- smoking.ever[ind1]
smk2 <- smoking.ever[ind2]
R.1 <- !(is.na(smk1))
R.2 <- !(is.na(smk2))

## Proportion of missing data.
sum(is.na(smk1)) / length(smk1)
sum(is.na(smk2)) / length(smk2)

## Complete-case analysis.
smok.cc <- summary(glm(smk2 ~ grs.std[ind2] + covariates[ind2, ], family = binomial))

## Inverse Probability Weighting.
smok.ipw <- ipw.logistic(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = smk2, R = R.2)$est

## Heckman's selection model.
smok.heck1 <- summary(selection(R.1 ~ grs.std[ind1] + covariates[ind1, ] + rct[ind1], as.logical(smk1) ~ grs.std[ind1] + covariates[ind1, ], method = "ml"))
smok.heck2 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs1.std[ind2], as.logical(smk2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
smok.heck3 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs2.std[ind2], as.logical(smk2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
smok.heck4 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs3.std[ind2], as.logical(smk2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))

## Tchetgen Tchetgen and Wirth.
smok.ttw.m1 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = smk1, hessian = TRUE)
smok.ttw.s1 <- sqrt(diag(solve(- smok.ttw.m1$hessian)))
smok.ttw.m2 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = smk2, hessian = TRUE)
smok.ttw.s2 <- sqrt(diag(solve(- smok.ttw.m2$hessian)))
smok.ttw.m3 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = smk2, hessian = TRUE)
smok.ttw.s3 <- sqrt(diag(solve(- smok.ttw.m3$hessian)))
smok.ttw.m4 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = smk2, hessian = TRUE)
smok.ttw.s4 <- sqrt(diag(solve(- smok.ttw.m4$hessian)))


## Assess standardized estimates.
smok.cc$coefficients
smok.ipw
smok.heck1
smok.heck2
smok.heck3
smok.heck4
cbind(smok.ttw.m1$par, smok.ttw.s1)
cbind(smok.ttw.m2$par, smok.ttw.s2)
cbind(smok.ttw.m3$par, smok.ttw.s3)
cbind(smok.ttw.m4$par, smok.ttw.s4)

## MR Estimates (using second-order approximation for standard errors).
c( smok.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(smok.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + smok.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
c( smok.ipw[2, 1] / bmi.ipw[2, 1], sqrt(smok.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + smok.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
c( smok.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(smok.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + smok.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
c( smok.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(smok.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + smok.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
c( smok.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(smok.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + smok.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
c( smok.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(smok.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + smok.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
c( smok.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(smok.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + smok.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
c( smok.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(smok.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + smok.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
c( smok.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(smok.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + smok.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
c( smok.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(smok.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + smok.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )

##################################################

##########   NUMBER OF CIGARETTES   ##########

## Number of cigarettes is a count variable, so we use 
## Poisson regression. We also consider linear regression 
## for comparison, since Heckman's method cannot work 
## with a Poisson model.

## Subset BMI as necessary.
smk1 <- smoking.freq[ind1]
smk2 <- smoking.freq[ind2]
R.1 <- !(is.na(smk1))
R.2 <- !(is.na(smk2))

## Proportion of missing data.
sum(is.na(smk1)) / length(smk1)
sum(is.na(smk2)) / length(smk2)



## First, run MR with a linear regression model for Y.

## Complete-case analysis.
ncig.l.cc <- summary(lm(smk2 ~ grs.std[ind2] + covariates[ind2, ]))

## Inverse Probability Weighting.
ncig.l.ipw <- ipw.linear(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = smk2, R = R.2)$est

## Heckman's selection model.
ncig.l.heck1 <- summary(selection(R.1 ~ grs.std[ind1] + covariates[ind1, ] + rct[ind1], smk1 ~ grs.std[ind1] + covariates[ind1, ], method = "2step"))
ncig.l.heck2 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs1.std[ind2], smk2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))
ncig.l.heck3 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs2.std[ind2], smk2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))
ncig.l.heck4 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs3.std[ind2], smk2 ~ grs.std[ind2] + covariates[ind2, ], method = "2step"))

## Tchetgen Tchetgen and Wirth.
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1])
ncig.l.ttw.m1 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = smk1, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.l.ttw.s1 <- sqrt(diag(solve(- ncig.l.ttw.m1$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2])
ncig.l.ttw.m2 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.l.ttw.s2 <- sqrt(diag(solve(- ncig.l.ttw.m2$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2])
ncig.l.ttw.m3 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.l.ttw.s3 <- sqrt(diag(solve(- ncig.l.ttw.m3$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2])
ncig.l.ttw.m4 <- optim(par = rep(0, 13), partial.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.l.ttw.s4 <- sqrt(diag(solve(- ncig.l.ttw.m4$hessian)))


## Assess standardized estimates.
ncig.l.cc$coefficients
ncig.l.ipw
ncig.l.heck1
ncig.l.heck2
ncig.l.heck3
ncig.l.heck4
cbind(ncig.l.ttw.m1$par, ncig.l.ttw.s1)
cbind(ncig.l.ttw.m2$par, ncig.l.ttw.s2)
cbind(ncig.l.ttw.m3$par, ncig.l.ttw.s3)
cbind(ncig.l.ttw.m4$par, ncig.l.ttw.s4)

## MR Estimates (using second-order approximation for standard errors).
c( ncig.l.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(ncig.l.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + ncig.l.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
c( ncig.l.ipw[2, 1] / bmi.ipw[2, 1], sqrt(ncig.l.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + ncig.l.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
c( ncig.l.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(ncig.l.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + ncig.l.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
c( ncig.l.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(ncig.l.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + ncig.l.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
c( ncig.l.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(ncig.l.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + ncig.l.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
c( ncig.l.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(ncig.l.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + ncig.l.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
c( ncig.l.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(ncig.l.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + ncig.l.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
c( ncig.l.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(ncig.l.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + ncig.l.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
c( ncig.l.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(ncig.l.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + ncig.l.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
c( ncig.l.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(ncig.l.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + ncig.l.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )



## Then run MR with a Poisson regression model for Y.

## Complete-case analysis.
ncig.cc <- summary(glm(smk2 ~ grs.std[ind2] + covariates[ind2, ], family = poisson))

## Inverse Probability Weighting.
ncig.ipw <- ipw.poisson(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = smk2, R = R.2)$est

## Heckman's method cannot be implemented for Poisson regression, 
## since the default R implementation does not support this model.

## Tchetgen Tchetgen and Wirth.
partial.opt21 <- optim(par = rep(0, 7), partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1])
ncig.ttw.m1 <- optim(par = rep(0, 12), partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = smk1, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.ttw.s1 <- sqrt(diag(solve(- ncig.ttw.m1$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2])
ncig.ttw.m2 <- optim(par = rep(0, 12), partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.ttw.s2 <- sqrt(diag(solve(- ncig.ttw.m2$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2])
ncig.ttw.m3 <- optim(par = rep(0, 12), partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.ttw.s3 <- sqrt(diag(solve(- ncig.ttw.m3$hessian)))
partial.opt21 <- optim(par = rep(0, 7), partial.poisson.lik1, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2])
ncig.ttw.m4 <- optim(par = rep(0, 12), partial.poisson.lik2, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = smk2, alpha.hat = partial.opt21$par, hessian = TRUE)
ncig.ttw.s4 <- sqrt(diag(solve(- ncig.ttw.m4$hessian)))


## Assess standardized estimates.
ncig.cc$coefficients
ncig.ipw
cbind(ncig.ttw.m1$par, ncig.ttw.s1)
cbind(ncig.ttw.m2$par, ncig.ttw.s2)
cbind(ncig.ttw.m3$par, ncig.ttw.s3)
cbind(ncig.ttw.m4$par, ncig.ttw.s4)

## MR Estimates (using second-order approximation for standard errors).
c( ncig.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(ncig.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + ncig.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
c( ncig.ipw[2, 1] / bmi.ipw[2, 1], sqrt(ncig.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + ncig.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
c( ncig.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(ncig.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + ncig.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
c( ncig.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(ncig.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + ncig.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
c( ncig.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(ncig.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + ncig.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
c( ncig.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(ncig.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + ncig.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )

##################################################

##########   SELF-HARM   ##########

## Self-harm is also binary, so we use logistic regression.

## Subset BMI as necessary.
har1 <- harm[ind1]
har2 <- harm[ind2]
R.1 <- !(is.na(har1))
R.2 <- !(is.na(har2))

## Proportion of missing data.
sum(is.na(har1)) / length(har1)
sum(is.na(har2)) / length(har2)

## Complete-case analysis.
harm.cc <- summary(glm(har2 ~ grs.std[ind2] + covariates[ind2, ], family = binomial))

## Inverse Probability Weighting.
harm.ipw <- ipw.logistic(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = har2, R = R.2)$est

## Heckman's selection model.
harm.heck1 <- summary(selection(R.1 ~ grs.std[ind1] + covariates[ind1, ] + rct[ind1], as.logical(har1) ~ grs.std[ind1] + covariates[ind1, ], method = "ml"))
harm.heck2 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs1.std[ind2], as.logical(har2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
harm.heck3 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs2.std[ind2], as.logical(har2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
harm.heck4 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs3.std[ind2], as.logical(har2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))

## Tchetgen Tchetgen and Wirth.
harm.ttw.m1 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = har1, hessian = TRUE)
harm.ttw.s1 <- sqrt(diag(solve(- harm.ttw.m1$hessian)))
harm.ttw.m2 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = har2, hessian = TRUE)
harm.ttw.s2 <- sqrt(diag(solve(- harm.ttw.m2$hessian)))
harm.ttw.m3 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = har2, hessian = TRUE)
harm.ttw.s3 <- sqrt(diag(solve(- harm.ttw.m3$hessian)))
harm.ttw.m4 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = har2, hessian = TRUE)
harm.ttw.s4 <- sqrt(diag(solve(- harm.ttw.m4$hessian)))


## Assess standardized estimates.
harm.cc$coefficients
harm.ipw
harm.heck1
harm.heck2
harm.heck3
harm.heck4
cbind(harm.ttw.m1$par, harm.ttw.s1)
cbind(harm.ttw.m2$par, harm.ttw.s2)
cbind(harm.ttw.m3$par, harm.ttw.s3)
cbind(harm.ttw.m4$par, harm.ttw.s4)

## MR Estimates (using second-order approximation for standard errors).
c( harm.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(harm.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + harm.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
c( harm.ipw[2, 1] / bmi.ipw[2, 1], sqrt(harm.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + harm.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
c( harm.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(harm.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + harm.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
c( harm.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(harm.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + harm.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
c( harm.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(harm.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + harm.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
c( harm.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(harm.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + harm.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
c( harm.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(harm.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + harm.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
c( harm.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(harm.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + harm.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
c( harm.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(harm.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + harm.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
c( harm.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(harm.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + harm.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )

##################################################

##########   DEPRESSION STATUS   ##########

## Depression is also a binary trait.

## Subset BMI as necessary.
dpr1 <- depression.all[ind1]
dpr2 <- depression.all[ind2]
R.1 <- !(is.na(dpr1))
R.2 <- !(is.na(dpr2))

## Proportion of missing data.
sum(is.na(dpr1)) / length(dpr1)
sum(is.na(dpr2)) / length(dpr2)

## Complete-case analysis.
deps.cc <- summary(glm(dpr2 ~ grs.std[ind2] + covariates[ind2, ], family = binomial))

## Inverse Probability Weighting.
deps.ipw <- ipw.logistic(X = cbind(grs.std[ind2], covariates[ind2, ]), Y = dpr2, R = R.2)$est

## Heckman's selection model.
deps.heck1 <- summary(selection(R.1 ~ grs.std[ind1] + covariates[ind1, ] + rct[ind1], as.logical(dpr1) ~ grs.std[ind1] + covariates[ind1, ], method = "ml"))
deps.heck2 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs1.std[ind2], as.logical(dpr2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
deps.heck3 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs2.std[ind2], as.logical(dpr2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))
deps.heck4 <- summary(selection(R.2 ~ grs.std[ind2] + covariates[ind2, ] + prs3.std[ind2], as.logical(dpr2) ~ grs.std[ind2] + covariates[ind2, ], method = "ml"))

## Tchetgen Tchetgen and Wirth.
deps.ttw.m1 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind1], covariates[ind1, ]), R = R.1, Z = rct[ind1], Y = dpr1, hessian = TRUE)
deps.ttw.s1 <- sqrt(diag(solve(- deps.ttw.m1$hessian)))
deps.ttw.m2 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs1.std[ind2], Y = dpr2, hessian = TRUE)
deps.ttw.s2 <- sqrt(diag(solve(- deps.ttw.m2$hessian)))
deps.ttw.m3 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs2.std[ind2], Y = dpr2, hessian = TRUE)
deps.ttw.s3 <- sqrt(diag(solve(- deps.ttw.m3$hessian)))
deps.ttw.m4 <- optim(par = rep(0, 19), full.logit.lik, method = "BFGS", control = list(fnscale = -1), X = cbind(grs.std[ind2], covariates[ind2, ]), R = R.2, Z = prs3.std[ind2], Y = dpr2, hessian = TRUE)
deps.ttw.s4 <- sqrt(diag(solve(- deps.ttw.m4$hessian)))


## Assess standardized estimates.
deps.cc$coefficients
deps.ipw
deps.heck1
deps.heck2
deps.heck3
deps.heck4
cbind(deps.ttw.m1$par, deps.ttw.s1)
cbind(deps.ttw.m2$par, deps.ttw.s2)
cbind(deps.ttw.m3$par, deps.ttw.s3)
cbind(deps.ttw.m4$par, deps.ttw.s4)

## MR Estimates (using second-order approximation for standard errors).
c( deps.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(deps.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + deps.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
c( deps.ipw[2, 1] / bmi.ipw[2, 1], sqrt(deps.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + deps.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
c( deps.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(deps.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + deps.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
c( deps.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(deps.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + deps.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
c( deps.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(deps.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + deps.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
c( deps.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(deps.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + deps.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
c( deps.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(deps.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + deps.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
c( deps.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(deps.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + deps.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
c( deps.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(deps.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + deps.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
c( deps.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(deps.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + deps.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )

##################################################

## As an additional analysis, we assess the strength 
## of each instrument for selection.

##################################################

##########  INSTRUMENT STRENGTH   ##########

## To quantify instrument strength, we use the R^2 and 
## F statistics with robust standard errors. These can
## be computed using the estimatr R package.


## First, for BMI.
bmi1 <- bmi[ind1]
bmi2 <- bmi[ind2]
R.1 <- !(is.na(bmi1))
R.2 <- !(is.na(bmi2))

## Proportion missing and prevalence.
sum(R.1 == 1) / length(R.1)
sum(R.2 == 1) / length(R.2)
mean(bmi1[R.1 == 1])
mean(bmi2[R.2 == 1])

## Compute R^2 statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$r.sq
summary(lm_robust(R.2 ~ prs1.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs2.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs3.std[ind2]))$r.sq

## Compute F statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$fstatistic
summary(lm_robust(R.2 ~ prs1.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs2.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs3.std[ind2]))$fstatistic



## For smoking status.
smk1 <- smoking.ever[ind1]
smk2 <- smoking.ever[ind2]
R.1 <- !(is.na(smk1))
R.2 <- !(is.na(smk2))

## Proportion missing and prevalence.
sum(R.1 == 1) / length(R.1)
sum(R.2 == 1) / length(R.2)
mean(smk1[R.1 == 1])
mean(smk2[R.2 == 1])

## Compute R^2 statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$r.sq
summary(lm_robust(R.2 ~ prs1.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs2.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs3.std[ind2]))$r.sq

## Compute F statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$fstatistic
summary(lm_robust(R.2 ~ prs1.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs2.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs3.std[ind2]))$fstatistic



## For number of cigarettes.
smk1 <- smoking.freq[ind1]
smk2 <- smoking.freq[ind2]
R.1 <- !(is.na(smk1))
R.2 <- !(is.na(smk2))

## Proportion missing and prevalence.
sum(R.1 == 1) / length(R.1)
sum(R.2 == 1) / length(R.2)
mean(smk1[R.1 == 1])
mean(smk2[R.2 == 1])

## Compute R^2 statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$r.sq
summary(lm_robust(R.2 ~ prs1.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs2.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs3.std[ind2]))$r.sq

## Compute F statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$fstatistic
summary(lm_robust(R.2 ~ prs1.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs2.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs3.std[ind2]))$fstatistic



## For self-harm.
har1 <- harm[ind1]
har2 <- harm[ind2]
R.1 <- !(is.na(har1))
R.2 <- !(is.na(har2))

## Proportion missing and prevalence.
sum(R.1 == 1) / length(R.1)
sum(R.2 == 1) / length(R.2)
mean(har1[R.1 == 1])
mean(har2[R.2 == 1])

## Compute R^2 statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$r.sq
summary(lm_robust(R.2 ~ prs1.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs2.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs3.std[ind2]))$r.sq

## Compute F statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$fstatistic
summary(lm_robust(R.2 ~ prs1.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs2.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs3.std[ind2]))$fstatistic



## For depression status.
dpr1 <- depression.all[ind1]
dpr2 <- depression.all[ind2]
R.1 <- !(is.na(dpr1))
R.2 <- !(is.na(dpr2))

## Proportion missing and prevalence.
sum(R.1 == 1) / length(R.1)
sum(R.2 == 1) / length(R.2)
mean(dpr1[R.1 == 1])
mean(dpr2[R.2 == 1])

## Compute R^2 statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$r.sq
summary(lm_robust(R.2 ~ prs1.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs2.std[ind2]))$r.sq
summary(lm_robust(R.2 ~ prs3.std[ind2]))$r.sq

## Compute F statistics.
summary(lm_robust(R.1 ~ rct[ind1]))$fstatistic
summary(lm_robust(R.2 ~ prs1.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs2.std[ind2]))$fstatistic
summary(lm_robust(R.2 ~ prs3.std[ind2]))$fstatistic

##################################################

## Now create LaTeX Tables and plots for the paper.

##########   RESULTS TABLES   ##########

## These are the actual Tables for the paper.

## First, for smoking and number of cigarettes (linear).
TableP1 <- matrix(NA, 10, 8)
colnames(TableP1) <- c("Smk Causal", "Smk StdErr", "Smk CI-lower", "Smk CI-upper", "N_cig Causal", "N_cig StdErr", "N_cig CI-lower", "N_cig CI-upper")
rownames(TableP1) <- c("CCA", "IPW", "Heck-RCT", "Heck-ALSPAC", "Heck-FFQ", "Heck-MHQ", "TTW-RCT", "TTW-ALSPAC", "TTW-FFQ", "TTW-MHQ")

## Smoking status.
TableP1[1, 1:2] <- c( smok.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(smok.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + smok.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
TableP1[2, 1:2] <- c( smok.ipw[2, 1] / bmi.ipw[2, 1], sqrt(smok.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + smok.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
TableP1[3, 1:2] <- c( smok.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(smok.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + smok.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
TableP1[4, 1:2] <- c( smok.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(smok.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + smok.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
TableP1[5, 1:2] <- c( smok.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(smok.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + smok.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
TableP1[6, 1:2] <- c( smok.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(smok.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + smok.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
TableP1[7, 1:2] <- c( smok.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(smok.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + smok.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
TableP1[8, 1:2] <- c( smok.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(smok.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + smok.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
TableP1[9, 1:2] <- c( smok.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(smok.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + smok.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
TableP1[10, 1:2] <- c( smok.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(smok.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + smok.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )
TableP1[, 3] <- TableP1[, 1] - qnorm(0.975) * TableP1[, 2]
TableP1[, 4] <- TableP1[, 1] + qnorm(0.975) * TableP1[, 2]

## Number of cigarettes.
TableP1[1, 5:6] <- c( ncig.l.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(ncig.l.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + ncig.l.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
TableP1[2, 5:6] <- c( ncig.l.ipw[2, 1] / bmi.ipw[2, 1], sqrt(ncig.l.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + ncig.l.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
TableP1[3, 5:6] <- c( ncig.l.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(ncig.l.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + ncig.l.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
TableP1[4, 5:6] <- c( ncig.l.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(ncig.l.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + ncig.l.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
TableP1[5, 5:6] <- c( ncig.l.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(ncig.l.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + ncig.l.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
TableP1[6, 5:6] <- c( ncig.l.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(ncig.l.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + ncig.l.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
TableP1[7, 5:6] <- c( ncig.l.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(ncig.l.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + ncig.l.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
TableP1[8, 5:6] <- c( ncig.l.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(ncig.l.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + ncig.l.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
TableP1[9, 5:6] <- c( ncig.l.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(ncig.l.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + ncig.l.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
TableP1[10, 5:6] <- c( ncig.l.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(ncig.l.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + ncig.l.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )
TableP1[, 7] <- TableP1[, 5] - qnorm(0.975) * TableP1[, 6]
TableP1[, 8] <- TableP1[, 5] + qnorm(0.975) * TableP1[, 6]


## Second, for self-harm and depression status.
TableP2 <- matrix(NA, 10, 8)
colnames(TableP2) <- c("Harm Causal", "Harm StdErr", "Harm CI-lower", "Harm CI-upper", "Depr Causal", "Depr StdErr", "Depr CI-lower", "Depr CI-upper")
rownames(TableP2) <- c("CCA", "IPW", "Heck-RCT", "Heck-ALSPAC", "Heck-FFQ", "Heck-MHQ", "TTW-RCT", "TTW-ALSPAC", "TTW-FFQ", "TTW-MHQ")

## Self-harm.
TableP2[1, 1:2] <- c( harm.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(harm.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + harm.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
TableP2[2, 1:2] <- c( harm.ipw[2, 1] / bmi.ipw[2, 1], sqrt(harm.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + harm.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
TableP2[3, 1:2] <- c( harm.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(harm.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + harm.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
TableP2[4, 1:2] <- c( harm.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(harm.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + harm.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
TableP2[5, 1:2] <- c( harm.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(harm.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + harm.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
TableP2[6, 1:2] <- c( harm.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(harm.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + harm.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
TableP2[7, 1:2] <- c( harm.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(harm.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + harm.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
TableP2[8, 1:2] <- c( harm.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(harm.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + harm.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
TableP2[9, 1:2] <- c( harm.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(harm.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + harm.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
TableP2[10, 1:2] <- c( harm.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(harm.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + harm.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )
TableP2[, 3] <- TableP2[, 1] - qnorm(0.975) * TableP2[, 2]
TableP2[, 4] <- TableP2[, 1] + qnorm(0.975) * TableP2[, 2]

## Depression status.
TableP2[1, 5:6] <- c( deps.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(deps.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + deps.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
TableP2[2, 5:6] <- c( deps.ipw[2, 1] / bmi.ipw[2, 1], sqrt(deps.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + deps.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
TableP2[3, 5:6] <- c( deps.heck1$estimate[9, 1] / bmi.heck1$estimate[9, 1], sqrt(deps.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^2 + deps.heck1$estimate[9, 1]^2 * bmi.heck1$estimate[9, 2]^2 / bmi.heck1$estimate[9, 1]^4) )
TableP2[4, 5:6] <- c( deps.heck2$estimate[9, 1] / bmi.heck2$estimate[9, 1], sqrt(deps.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^2 + deps.heck2$estimate[9, 1]^2 * bmi.heck2$estimate[9, 2]^2 / bmi.heck2$estimate[9, 1]^4) )
TableP2[5, 5:6] <- c( deps.heck3$estimate[9, 1] / bmi.heck3$estimate[9, 1], sqrt(deps.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^2 + deps.heck3$estimate[9, 1]^2 * bmi.heck3$estimate[9, 2]^2 / bmi.heck3$estimate[9, 1]^4) )
TableP2[6, 5:6] <- c( deps.heck4$estimate[9, 1] / bmi.heck4$estimate[9, 1], sqrt(deps.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^2 + deps.heck4$estimate[9, 1]^2 * bmi.heck4$estimate[9, 2]^2 / bmi.heck4$estimate[9, 1]^4) )
TableP2[7, 5:6] <- c( deps.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(deps.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + deps.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
TableP2[8, 5:6] <- c( deps.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(deps.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + deps.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
TableP2[9, 5:6] <- c( deps.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(deps.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + deps.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
TableP2[10, 5:6] <- c( deps.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(deps.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + deps.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )
TableP2[, 7] <- TableP2[, 5] - qnorm(0.975) * TableP2[, 6]
TableP2[, 8] <- TableP2[, 5] + qnorm(0.975) * TableP2[, 6]


## Third, for number of cigarettes (Poisson).
TableP3 <- matrix(NA, 6, 4)
colnames(TableP3) <- c("N_cig Causal", "N_cig StdErr", "N_cig CI-lower", "N_cig CI-upper")
rownames(TableP3) <- c("CCA", "IPW", "TTW-RCT", "TTW-ALSPAC", "TTW-FFQ", "TTW-MHQ")

## Number of cigarettes.
TableP3[1, 1:2] <- c( ncig.cc$coef[2, 1] / bmi.cc$coef[2, 1], sqrt(ncig.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^2 + ncig.cc$coef[2, 1]^2 * bmi.cc$coef[2, 2]^2 / bmi.cc$coef[2, 1]^4) )
TableP3[2, 1:2] <- c( ncig.ipw[2, 1] / bmi.ipw[2, 1], sqrt(ncig.ipw[2, 2]^2 / bmi.ipw[2, 1]^2 + ncig.ipw[2, 1]^2 * bmi.ipw[2, 2]^2 / bmi.ipw[2, 1]^4) )
TableP3[3, 1:2] <- c( ncig.ttw.m1$par[2] / bmi.ttw.m1$par[2], sqrt(ncig.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^2 + ncig.ttw.m1$par[2]^2 * bmi.ttw.s1[2]^2 / bmi.ttw.m1$par[2]^4) )
TableP3[4, 1:2] <- c( ncig.ttw.m2$par[2] / bmi.ttw.m2$par[2], sqrt(ncig.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^2 + ncig.ttw.m2$par[2]^2 * bmi.ttw.s2[2]^2 / bmi.ttw.m2$par[2]^4) )
TableP3[5, 1:2] <- c( ncig.ttw.m3$par[2] / bmi.ttw.m3$par[2], sqrt(ncig.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^2 + ncig.ttw.m3$par[2]^2 * bmi.ttw.s3[2]^2 / bmi.ttw.m3$par[2]^4) )
TableP3[6, 1:2] <- c( ncig.ttw.m4$par[2] / bmi.ttw.m4$par[2], sqrt(ncig.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^2 + ncig.ttw.m4$par[2]^2 * bmi.ttw.s4[2]^2 / bmi.ttw.m4$par[2]^4) )
TableP3[, 3] <- TableP3[, 1] - qnorm(0.975) * TableP3[, 2]
TableP3[, 4] <- TableP3[, 1] + qnorm(0.975) * TableP3[, 2]


## Report them in format that can be transferred to LaTeX.
xtable(TableP1, digits = 3)
xtable(TableP2, digits = 3)
xtable(TableP3, digits = 3)

###################################################

##########   FOREST PLOTS   ##########

## It may be prettier to report results in terms of 
## forest plots, rather than tables.

## So here we go, for smoking and number of cigarettes.

## Start plotting.
pdf(file = "ALSPAC_bmi_smk.pdf", width = 14, height = 6, pointsize = 16)

## Set up main parameters.
par(mfrow = c(1, 2))
par(mar = c(5, 8, 3, 0))

## Set up colors.
colors <- rev(c("black", "black", NA, "red", "blue", "darkgreen", "purple", NA, "red", "blue", "darkgreen", "purple"))

## Plot 1 arguments.
plot1_means <- rev(c(TableP1[1:2, 1], NA, TableP1[3:6, 1], NA, TableP1[7:10, 1]))
plot1_ci1 <- rev(c(TableP1[1:2, 3], NA, TableP1[3:6, 3], NA, TableP1[7:10, 3]))
plot1_ci2 <- rev(c(TableP1[1:2, 4], NA, TableP1[3:6, 4], NA, TableP1[7:10, 4]))
plot1_traits <- rev(c("CCA", "IPW", "Heck (RCT)", "Heck (ALSPAC)", "Heck (FFQ)", "Heck (MHQ)", "TTW (RCT)", "TTW (ALSPAC)", "TTW (FFQ)", "TTW (MHQ)"))

## Do plot 1.
plot(x = plot1_means, y = 1:12, type = "p", axes = FALSE, main = "Smoking (Ever vs Never)", xlab = "Causal Effect", ylab = "", xlim = c(-0.4, 0.4), ylim = c(0.5, 12), pch = 19, cex.lab = 1.1, cex.main = 1.1, col = colors)
for (i in 1:12) {
  lines(c(plot1_ci1[i], plot1_ci2[i]), c(i, i), col = colors[i])
  lines(c(plot1_ci1[i], plot1_ci1[i]), c(i - 0.1, i + 0.1), col = colors[i])
  lines(c(plot1_ci2[i], plot1_ci2[i]), c(i - 0.1, i + 0.1), col = colors[i])
}
axis(side = 2, at = c(1:4, 6:9, 11:12), labels = plot1_traits, las = 1, cex.axis = 0.9)
axis(side = 1, at = c(-0.4, -0.2, 0, 0.2, 0.4), las = 1)
abline(v = c(-0.4, -0.2, 0, 0.2, 0.4), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")

## Update main parameters.
par(mar = c(5, 4, 3, 2))

## Plot 2 arguments.
plot2_means <- rev(c(TableP1[1:2, 5], NA, TableP1[3:6, 5], NA, TableP1[7:10, 5]))
plot2_ci1 <- rev(c(TableP1[1:2, 7], NA, TableP1[3:6, 7], NA, TableP1[7:10, 7]))
plot2_ci2 <- rev(c(TableP1[1:2, 8], NA, TableP1[3:6, 8], NA, TableP1[7:10, 8]))
#plot2_traits <- rev(c("CCA", "IPW", "Heck (RCT)", "Heck (ALSPAC)", "Heck (FFQ)", "Heck (MHQ)", "TTW (RCT)", "TTW (ALSPAC)", "TTW (FFQ)", "TTW (MHQ)"))

## Do plot 2.
plot(x = plot2_means, y = 1:12, type = "p", axes = FALSE, main = "Number of Cigarettes per Day", xlab = "Causal Effect", ylab = "", xlim = c(-2, 4), ylim = c(0.5, 12), pch = 19, cex.lab = 1.1, cex.main = 1.1, col = colors)
for (i in 1:12) {
  lines(c(plot2_ci1[i], plot2_ci2[i]), c(i, i), col = colors[i])
  lines(c(plot2_ci1[i], plot2_ci1[i]), c(i - 0.1, i + 0.1), col = colors[i])
  lines(c(plot2_ci2[i], plot2_ci2[i]), c(i - 0.1, i + 0.1), col = colors[i])
}
#axis(side = 2, at = c(1:4, 6:9, 11:12), labels = plot2_traits, las = 1, cex.axis = 0.9)
axis(side = 1, at = c(-2, 0, 2, 4), las = 1)
abline(v = c(-2, 0, 2, 4), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")

## Goodbye.
dev.off()



## Second plot, for depression and self-harm.

## Start plotting.
pdf(file = "ALSPAC_bmi_mental.pdf", width = 14, height = 6, pointsize = 16)

## Set up main parameters.
par(mfrow = c(1, 2))
par(mar = c(5, 8, 3, 0))

## Set up colors.
colors <- rev(c("black", "black", NA, "red", "blue", "darkgreen", "purple", NA, "red", "blue", "darkgreen", "purple"))

## Plot 1 arguments.
plot1_means <- rev(c(TableP2[1:2, 1], NA, TableP2[3:6, 1], NA, TableP2[7:10, 1]))
plot1_ci1 <- rev(c(TableP2[1:2, 3], NA, TableP2[3:6, 3], NA, TableP2[7:10, 3]))
plot1_ci2 <- rev(c(TableP2[1:2, 4], NA, TableP2[3:6, 4], NA, TableP2[7:10, 4]))
plot1_traits <- rev(c("CCA", "IPW", "Heck (RCT)", "Heck (ALSPAC)", "Heck (FFQ)", "Heck (MHQ)", "TTW (RCT)", "TTW (ALSPAC)", "TTW (FFQ)", "TTW (MHQ)"))

## Do plot 1.
plot(x = plot1_means, y = 1:12, type = "p", axes = FALSE, main = "Self-Harm Status", xlab = "Causal Effect", ylab = "", xlim = c(-0.5, 0.5), ylim = c(0.5, 12), pch = 19, cex.lab = 1.1, cex.main = 1.1, col = colors)
for (i in 1:12) {
  lines(c(plot1_ci1[i], plot1_ci2[i]), c(i, i), col = colors[i])
  lines(c(plot1_ci1[i], plot1_ci1[i]), c(i - 0.1, i + 0.1), col = colors[i])
  lines(c(plot1_ci2[i], plot1_ci2[i]), c(i - 0.1, i + 0.1), col = colors[i])
}
axis(side = 2, at = c(1:4, 6:9, 11:12), labels = plot1_traits, las = 1, cex.axis = 0.9)
axis(side = 1, at = c(-0.5, -0.25, 0, 0.25, 0.5), las = 1)
abline(v = c(-0.5, -0.25, 0, 0.25, 0.5), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")

## Update main parameters.
par(mar = c(5, 4, 3, 2))

## Plot 2 arguments.
plot2_means <- rev(c(TableP2[1:2, 5], NA, TableP2[3:6, 5], NA, TableP2[7:10, 5]))
plot2_ci1 <- rev(c(TableP2[1:2, 7], NA, TableP2[3:6, 7], NA, TableP2[7:10, 7]))
plot2_ci2 <- rev(c(TableP2[1:2, 8], NA, TableP2[3:6, 8], NA, TableP2[7:10, 8]))
#plot2_traits <- rev(c("CCA", "IPW", "Heck (RCT)", "Heck (ALSPAC)", "Heck (FFQ)", "Heck (MHQ)", "TTW (RCT)", "TTW (ALSPAC)", "TTW (FFQ)", "TTW (MHQ)"))

## Do plot 2.
plot(x = plot2_means, y = 1:12, type = "p", axes = FALSE, main = "Depression Status", xlab = "Causal Effect", ylab = "", xlim = c(-1, 1), ylim = c(0.5, 12), pch = 19, cex.lab = 1.1, cex.main = 1.1, col = colors)
for (i in 1:12) {
  lines(c(plot2_ci1[i], plot2_ci2[i]), c(i, i), col = colors[i])
  lines(c(plot2_ci1[i], plot2_ci1[i]), c(i - 0.1, i + 0.1), col = colors[i])
  lines(c(plot2_ci2[i], plot2_ci2[i]), c(i - 0.1, i + 0.1), col = colors[i])
}
#axis(side = 2, at = c(1:4, 6:9, 11:12), labels = plot2_traits, las = 1, cex.axis = 0.9)
axis(side = 1, at = c(-1, -0.5, 0, 0.5, 1), las = 1)
abline(v = c(-1, -0.5, 0, 0.5, 1), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")

## Goodbye.
dev.off()



## Third plot, for number of cigarettes using Poisson regression.

## Start plotting.
pdf(file = "ALSPAC_bmi_cig2.pdf", width = 7, height = 4, pointsize = 16)

## Set up main parameters.
par(mfrow = c(1, 1))   ## This will be a single plot.
par(mar = c(5, 8, 3, 1))

## Set up colors.
colors <- rev(c("black", "black", NA, "red", "blue", "darkgreen", "purple"))

## Plot 1 arguments (no Heckman here).
plot1_means <- rev(c(TableP3[1:2, 1], NA, TableP3[3:6, 1]))
plot1_ci1 <- rev(c(TableP3[1:2, 3], NA, TableP3[3:6, 3]))
plot1_ci2 <- rev(c(TableP3[1:2, 4], NA, TableP3[3:6, 4]))
plot1_traits <- rev(c("CCA", "IPW", "TTW (RCT)", "TTW (ALSPAC)", "TTW (FFQ)", "TTW (MHQ)"))

## Do plot 1.
plot(x = plot1_means, y = 1:7, type = "p", axes = FALSE, main = "Number of Cigarettes per Day", xlab = "Causal Effect", ylab = "", xlim = c(-1, 3), ylim = c(0.5, 7), pch = 20, cex.lab = 1.1, cex.main = 1.1, col = colors)
for (i in 1:7) {
  lines(c(plot1_ci1[i], plot1_ci2[i]), c(i, i), col = colors[i])
  lines(c(plot1_ci1[i], plot1_ci1[i]), c(i - 0.15, i + 0.15), col = colors[i])
  lines(c(plot1_ci2[i], plot1_ci2[i]), c(i - 0.15, i + 0.15), col = colors[i])
}
axis(side = 2, at = c(1:4, 6:7), labels = plot1_traits, las = 1, cex.axis = 0.9)
axis(side = 1, at = c(-1, 0, 1, 2, 3), las = 1)
abline(v = c(-1, 0, 1, 2, 3), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")

## No Plot 2 here.

## Goodbye.
dev.off()

###################################################




##################################################
