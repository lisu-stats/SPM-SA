# SPM-SA
Code for performing  transparent sensitivity analysis for shared parameter models, based on the manuscript 'A Sensitivity Analysis Approach for Informative Dropout using Shared Parameter Models' by Li Su, Qiuju Li, Jessica Barrett and Mike Daniels (2018).

1) SPM_model.R: WinBUGS code for fitting the shared parameter model to the HERS data
2) postsample100.Rdata: R data file for saved 100 posterior samples from the fitted SPM to the HERS data
3) EV_skewednormaldist.R: R functions to calcuate the mean and variance of a multi-variate skew-normal distribution; to be called when running G-compuation in sensitivity analyses
4) G-computation.R: R function to run G-compuation based on one posterior sample from the fitted SPM; R code to run the G-compuation in parallel
5) summary-G-computation.R: R code to summarize the G-compuation results in the HERS analysis
