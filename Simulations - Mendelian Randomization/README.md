
# Simulations - MR

This folder contains the R files from the simulation study to assess the performance of "instruments for selection" methods in adjusting for bias due to missing data in Mendelian Randomization. In total, we considered 28 simulation scenarios, which are briefly described below:

 - Scenarios 1-6 contain simulations for one-sample MR analyses with a single genetic instrument for inference (possibly a polygenic score). 
 - Scenarios 7-12 contain similar simulations as scenarios 1-6, but for two-sample MR analyses instead. 
 - Scenarios 13-18 contain simulations for one-sample MR analyses with multiple instruments for inference.
 - Scenarios 19-24 contain similar simulations as scenarios 13-18, but for two-sample MR analyses instead. 

Each of these categories contains six simulation scenarios, corresponding to either a null or a positive causal effect, and missingness affected by either the exposure or the risk factor or both.

 - Finally, scenarios 25-28 contain a brief exploration of the performance of IVsel methods in simulations with a misspecified instrument-exposure model and/or population differences.

We refer to our manuscript for more details.

The code used to run each simulation is included in the files "MR_Sim##.R", while simulation results are reported in files "MR##_results.RData". The implementation of the IVsel methods for these simulations was made using the R functions contained in the files "IVsel_functions".

Finally, the file "Summary.R" explores the simulation results from the various scenarios and describes how they were combined to form the Tables included in our manuscript and its supplement.

