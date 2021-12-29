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
- To run the algorithm, simply input the time-to-event outcome, status, covariate matrix, and binary treatment to the R function TBAFTcure and specify the number of iterations. A toy example is given below and further details can be find in simulation.R.
```
n = 1000
p = 8
checktime=1
x = matrix(rnorm(n*p), nrow=n) 
x[,5] = rbinom(n,1,0.5)
x[,4] = sample(1:3,n,replace = TRUE, prob = c(0.2,0.4,0.4))
g = 2*(x[,4]==1) - (x[,4]==2) - 4*(x[,4] == 3)
x[,8] = runif(n,-1,0)
q = -6 + g + 6*abs(x[,3]-1)
  
# generate treatment variable
pi = 0.8*pnorm(3*q/sd(q) - 0.5*x[,1]) + 0.05 + runif(n)/10
a = rbinom(n,1,pi)
# the cured fraction
cured_p = 0.8*pnorm(- 0.5*x[,7] + x[,4]*a) + 0.05 - 0.1*x[,8] 
true_label = rbinom(n,1,cured_p)
tau = 1 + .25*x[,2]*x[,5] + 0.5*x[,2]*(1-x[,5])
mu = (0.2*q + tau*a)

sigma = 0.5
## draw the response variable with additive error
y.true = y = exp(mu + sigma*rnorm(n))
y.true[which(true_label==0)] = y[which(true_label==0)] = 1e10
c = rexp(n, 0.5)
delta = (y.true <= c)
delta[true_label==0]=FALSE
y[delta == 0] = c[delta == 0]
#pihat = glm(a ~ x, family = "binomial"(link = "probit"))$fitted.values
bartfit <- BART::pbart(x, a)
pihat <- as.matrix(bartfit$prob.train.mean)
cat('Number of uncured: ', sum(true_label), ' .\n')
cat('Number of failed: ', sum(delta), ' .\n')
fit <- TBAFTcure(y, delta, a, x, x, x_cure = cbind(x,a), pihat)
```
- Outputs:
 * posterior draws of (U)CATE: fit$tau
 * posterior draws of group label: fit$cured_label
 * variable usage count for control/modifer/cured tree: fit$varcnt_con, fit$varcnt_mod, fit$varcnt_cure
