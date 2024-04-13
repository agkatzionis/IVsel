
# Simulations - MR

This folder contains the R files from the simulation study to assess the performance of "instruments for selection" methods in adjusting for bias due to missing data in Mendelian Randomization. In total, we considered 28 simulation scenarios, which are briefly described below:

 - Scenarios 1-6 contain simulations for one-sample MR analyses with a single genetic instrument for inference (possibly a polygenic score). 
 - Scenarios 7-12 contain similar simulations as scenarios 1-6, but for two-sample MR analyses instead. 
 - Scenarios 13-18 contain simulations for one-sample MR analyses with multiple instruments for inference.
 - Scenarios 19-24 contain similar simulations as scenarios 13-18, but for two-sample MR analyses instead. 

Each of these categories contains six simulation scenarios, corresponding to either a null or a positive causal effect, and missingness affected by either the exposure or the risk factor or both.

 - Finally, scenarios 25-28 contain a brief exploration of the performance of IVsel methods in simulations with a misspecified instrument-exposure model and/or population differences.

We refer to our manuscript for more details.

The code used to run each simulation is included in the files "MR_Sim##.R", while simulation results are reported in files "MR##_results.RData". The implementation of the IVsel methods for these simulations was made using the R functions contained in the files "IVsel_functions". The file "Summary.R" explores the simulation results from the various scenarios and describes how they were combined to form the Tables and Figures included in our manuscript and its supplement.

Finally, after the manuscript had been submitted for publication, some additional simulation scenarios were added in response to reviewers' comments. For brevity, these scenarios were only simulated under a one-cample MR setting with a single instrument for inference, in an MR analysis where both the exposure and the outcome affected selection into the study. The following scenarios were considered:

 - File "MrSimR1_1_XY.R" contains MR simulations with a varying proportion of missing values (6 scenarios). Files "MrSimR1_1_XY_#_results.RData" contain the results.
 - File "MrSimR1_2_XY.R" contains MR simulations with an added effect of the instrument for selection on the MR exposure, of varyinng strength (5 scenarios). Files "MrSimR1_2_XY_#_results.RData" contain the results.
 - File "MrSimR1_3_XY.R" contains MR simulations with an added effect of the instrument for selection on the MR outcome, of varyinng strength (5 scenarios). Files "MrSimR1_3_XY_#_results.RData" contain the results.
 - File "MrSimR1_4_XY.R" contains MR simulations with different values for the MR causal effect (6 scenarios). Files "MrSimR1_4_XY_#_results.RData" contain the results.
 - File "MrSimR1_5_XY.R" contains MR simulations with different sample sizes (6 scenarios). Files "MrSimR1_5_XY_#_results.RData" contain the results.

The file "Summary.R" was updated to include simulation results from these 28 additional scenarios.
