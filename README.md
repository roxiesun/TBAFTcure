# Code for TBAFTcure
Authors: Rongqian Sun and Xinyuan Song

Paper link: tbc

Contact: <sunrq@link.cuhk.edu.hk>

This directory provides implementation of the tree-based Bayesian AFT cure model for heterogenous treatment effect estimation with simulated dataset.

## Implementation details
```
packages <- c('Rcpp', 'RcppArmadillo','survival', 'survminer')
lapply(packages, require, character.only = TRUE)
source('TBAFTcure.cpp')
source('TBAFTcure.R')
```
