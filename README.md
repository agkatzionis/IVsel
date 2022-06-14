# Instruments for Selection

This directory contains R code in support of the paper "Using Instruments for Selection to Adjust for Selection Bias in Mendelian Randomization". It includes implementations of inverse probability weighting, Heckman's sample selection model and the selection bias adjustment method of Tchetgen Tchetgen and Wirth (2017). These methods are implemented for regression analyses (linear, logistic and Poisson regression) as well as for Mendelian randomization studies. R code to compute the selection bounds of Marden et al. (2018) is also provided.

In addition to the R functions that can be used to implement these methods, this repository also includes the R code used to conduct the simulation study included in the original paper, as well as some of the code used to run the real-data application using ALSPAC data (the data are not available here for privacy reasons).
