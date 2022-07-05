
# Simulations - Regression Analysis

This folder contains the R files from the simulation study to assess the performance of "instruments for selection" methods in adjusting for bias due to missing data in regression-type simulations. In total, this folder includes 50 simulation scenarios. Here is a brief summary of each scenario:

The results of the first 15 scenarios are reported in the main part of our paper:

 - Scenario 1 is the baseline scenario for our simulation study.
 - Scenarios 2-10 explore a range of simulation settings and are the scenarios reported in Table 1 of our paper. They include a null exposure outcome association (Scen 2), a binary instrument (Scen 3), a binary exposure (Scen 4), a binary outcome (Scen 5), a Poisson-distributed outcome (Scen 6), a selection process affected by only the outcome (Scen 7) or only the exposure (Scen 8), and an instrument exerting direct effects on the exposure (Scen 9) or the outcome (Scen 10).
 - Scenarios 11-15 explore th eperformance of IVsel methods when we vary the instrument strength of the instrument for selection.
 
 The results of the remaining 35 scenarios are reported in the supplement of our paper:
 
 - Scenarios S1-S8 replicate the simulations of Scenarios 3-10 but for a null exposure-outcome association (instead of a positive one).
 - Scenarios S9-S14 explore the performance of IVsel methods when varying the proportion of individuals with fully observed data in our sample.
 - Scenarios S15-S19 explore the performance of IVsel methods when varying the strength of the exposure-selection effect.
 - Scenario S20 considers the case where the instrument for selection is in fact not associated with selection, thus violating the relevance assumption.
 - Scenarios S21-S25 explore the performance of IVsel methods when varying the strength of the outcome-selection effect.
 - Scenarios S26-S27 contain simulations with a weak binary instrument.
 - Scenarios 28-29 explore the performance of IVsel methods when the outcome does not cause selection, but instead the two are confounded.
 - Scenarios 30-35 consider violations of the normality assumption about the outcome.

We refer to our manuscript for more details.

The code used to run each simulation is included in the files "Reg_Sim##.R", while simulation results are reported in files "Reg###_results.RData". The implementation of the IVsel methods for these simulations was made using the R functions contained in the files "IVsel_functions".

Finally, the file "Summary.R" explores the simulation results from the various scenarios and describes how they were combined to form the Tables included in our manuscript and its supplement.
