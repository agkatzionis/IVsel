# Instruments for Selection

This directory contains R code in support of the paper "Using Instruments for Selection to Adjust for Selection Bias in Mendelian Randomization". It includes implementations of inverse probability weighting and the selection bias adjustment method of [Tchetgen Tchetgen and Wirth (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12670). These methods are implemented for regression analyses (mainly linear, logistic and Poisson regression), but can easily be adapted for Mendelian randomization studies. A way to compute the selection bounds of [Marden et al. (2018)](https://journals.lww.com/epidem/Fulltext/2018/05000/Implementation_of_Instrumental_Variable_Bounds_for.7.aspx) is also provided.

The R code to implement these functions is given in the file "IVsel_Code.R", along with some instructions on how to use them. The functions are also included in the R workspace "IVsel_RData".

In addition, the subfolders of this repository contain the R code used to conduct the simulation study included in the paper, as well as some of the code used to run the real-data stuy on the effects of high BMI on smoking traits in early adolescence using [ALSPAC](http://www.bristol.ac.uk/alspac/) data (the raw data are not available here for privacy reasons).
