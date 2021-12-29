# Code for TBAFTcure
Authors: Rongqian Sun and Xinyuan Song

Paper link: tbc

Contact: <sunrq@link.cuhk.edu.hk>

This directory provides implementation of the tree-based Bayesian AFT cure model for heterogenous treatment effect estimation with simulated dataset.

## Implementation details
```
packages <- c('Rcpp', 'RcppArmadillo','survival', 'survminer')
lapply(packages, require, character.only = TRUE)
source('TBAFTcure.cpp')tr
## or
# system("R CMD INSTALL TBAFTcure_0.1.0.tar.gz")
# packages <- c('Rcpp', 'RcppArmadillo','survival', 'survminer','TBAFTcure')
# lapply(packages, require, character.only = TRUE)
```
- To run the algorithm, simply input the time-to-event outcome, status, covariate matrix, and binary treatment to the R function TBAFTcure and specify the number of iterations. 
