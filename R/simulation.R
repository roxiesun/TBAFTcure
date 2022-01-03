packages <- c('Rcpp', 'RcppArmadillo', 'mvtnorm', 'survival', 'AFTrees') #, 'TBAFTcure')
lapply(packages, require, character.only = TRUE)
sourceCpp('TBAFTcure.cpp')
source('TBAFTcure.R')
#data generating process
n = 1000
p = 8
checktime=1
nrep = 1
flag = 1
startseed = 1001

res_95 = res_99 = res_obj = res_objsp = array(NA,dim = c(nrep,5))
colnames(res_95) = colnames(res_99) = colnames(res_obj) = colnames(res_objsp) = c('bias.UCATE','RMS.UCATE','avg.cov.UCATE','bias.sigma', 'AVG.len')

SATE_95 = SATE_99 = SATE_obj = SATE_objsp = array(NA,dim = c(nrep,4))
colnames(SATE_95) = colnames(SATE_99) = colnames(SATE_obj) = colnames(SATE_objsp) = c('bias.USATE','RMS_comp.USATE','cov.USATE', 'AVG.len')

SurvPr_95 = SurvPr_99 = SurvPr_obj = SurvPr_objsp = array(NA, dim = c(nrep,4))
colnames(SurvPr_95) = colnames(SurvPr_99) = colnames(SurvPr_obj) = colnames(SurvPr_objsp) = c('bias.CASP', 'RMS.CASP', 'avg.cov.CASP', 'AVG.len')

time.begin <- proc.time()
for (irep in 1:nrep) {
  set.seed(startseed + irep - 1)
  
  x = matrix(rnorm(n*p), nrow=n) 
  x[,5] = rbinom(n,1,0.5)
  x[,4] = sample(1:3,n,replace = TRUE, prob = c(0.2,0.4,0.4))
  g = 2*(x[,4]==1) - (x[,4]==2) - 4*(x[,4] == 3)
  x[,8] = runif(n,-1,0)
  
  # create targeted selection
  q = -6 + g + 6*abs(x[,3]-1)
  
  # generate treatment variable
  #pi = pnorm(q)
  pi = 0.8*pnorm(3*q/sd(q) - 0.5*x[,1]) + 0.05 + runif(n)/10
  a = rbinom(n,1,pi)
  
  # the cured fraction
  u = runif(n)/10
  cured_p = 0.8*pnorm(- 0.5*x[,7] + x[,4]*a) + 0.05 - 0.1*x[,8] 
  p_1 = 0.8*pnorm(- 0.5*x[,7] + x[,4]) + 0.05 - 0.1*x[,8] 
  p_0 = 0.8*pnorm(- 0.5*x[,7]) + 0.05 -0.1*x[,8]
  
  true_label = rbinom(n,1,cured_p)
  
  tau = 1 + .25*x[,2]*x[,5] + 0.5*x[,2]*(1-x[,5])
  mu = (0.2*q + tau*a)
  # set the noise level relative to the expected mean function of Y
  sigma = 0.5
  ## draw the response variable with additive error
  ## setting1: normal error
  y.true = y = exp(mu + sigma*rnorm(n))
  s_1 = pnorm((log(checktime) - 0.2*q - tau)/sigma, lower.tail = FALSE)
  s_0 = pnorm((log(checktime) - 0.2*q)/sigma, lower.tail = FALSE)
  y.true[which(true_label==0)] = y[which(true_label==0)] = 1e10
  
  true_survp = p_0 - p_1 + p_1*s_1 - p_0*s_0
  c = rexp(n, 0.5)
  delta = (y.true <= c)
  delta[true_label==0]=FALSE
  y[delta == 0] = c[delta == 0]
  pihat = glm(a ~ x, family = "binomial"(link = "probit"))$fitted.values #npihat = BART::pbart(x,a)
  
  
  
  AFTbcf_fit1 = TBAFTcure(y, delta, a, x, x_cure = cbind(x,a), pihat = pihat,
                               nburn=1000, nsim=2000, nthin = 2, ntree_control = 200,
                               ntree_moderate = 50,ntree_cure = 200, sigq = 0.95)


  AFTbcf_fit2 = TBAFTcure(y, delta, a, x, x_cure = cbind(x,a), pihat = pihat,
                               nburn=1000, nsim=2000, nthin = 2, ntree_control = 200,
                               ntree_moderate = 50,ntree_cure = 200,sigq = 0.99)
  
  
  obj <- AFTrees::IndivAFT(x.train=cbind(x,pihat),y.train=y,status=delta,Trt=a, 
                  ndpost=4000,nskip=1000,printevery=100, sigquant=.5, keepevery=2)
  
  obj_sp <- AFTrees::IndivAFT(x.train=cbind(x,pihat),y.train=y,status=delta,Trt=a,
                     nonparametric=FALSE, ndpost=4000,nskip=1000,printevery=100,
                     sigquant=.9, keepevery=2)

  
  res_95[irep,1] = mean(colMeans(AFTbcf_fit1$tau)[true_label==1] - tau[true_label==1])
  res_95[irep,2] = sqrt(mean(((colMeans(AFTbcf_fit1$tau) - tau)[true_label==1])^2))
  l.ci.sp <- apply(AFTbcf_fit1$tau[,true_label==1], 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(AFTbcf_fit1$tau[,true_label==1], 2, function(x) quantile(x, prob=.975))
  res_95[irep,3] = mean((tau[true_label==1] > l.ci.sp) & (tau[true_label==1] <= u.ci.sp))
  res_95[irep,4] = mean(AFTbcf_fit1$sigma) - sigma
  res_95[irep,5] = mean(u.ci.sp - l.ci.sp)

  res_99[irep,1] = mean(colMeans(AFTbcf_fit2$tau)[true_label==1] - tau[true_label==1])
  res_99[irep,2] = sqrt(mean(((colMeans(AFTbcf_fit2$tau) - tau)[true_label==1])^2))
  l.ci.sp <- apply(AFTbcf_fit2$tau[,true_label==1], 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(AFTbcf_fit2$tau[,true_label==1], 2, function(x) quantile(x, prob=.975))
  res_99[irep,3] = mean((tau[true_label==1] > l.ci.sp) & (tau[true_label==1] <= u.ci.sp))
  res_99[irep,4] = mean(AFTbcf_fit2$sigma) - sigma
  res_99[irep,5] = mean(u.ci.sp - l.ci.sp)
  
  res_obj[irep,1] = mean(colMeans(obj$Theta)[true_label==1] - tau[true_label==1])
  res_obj[irep,2] = sqrt(mean(((colMeans(obj$Theta) - tau)[true_label==1])^2))
  l.ci.sp <- apply(obj$Theta[,true_label==1], 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(obj$Theta[,true_label==1], 2, function(x) quantile(x, prob=.975))
  res_obj[irep,3] = mean((tau[true_label==1] > l.ci.sp) & (tau[true_label==1] <= u.ci.sp))
  res_obj[irep,4] = mean(obj$sigma) - sigma
  res_obj[irep,5] = mean(u.ci.sp - l.ci.sp)
  
  res_objsp[irep,1] = mean(colMeans(obj_sp$Theta)[true_label==1] - tau[true_label==1])
  res_objsp[irep,2] = sqrt(mean(((colMeans(obj_sp$Theta) - tau)[true_label==1])^2))
  l.ci.sp <- apply(obj_sp$Theta[,true_label==1], 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(obj_sp$Theta[,true_label==1], 2, function(x) quantile(x, prob=.975))
  res_objsp[irep,3] = mean((tau[true_label==1] > l.ci.sp) & (tau[true_label==1] <= u.ci.sp))
  res_objsp[irep,4] = mean(obj_sp$sigma) - sigma
  res_objsp[irep,5] = mean(u.ci.sp - l.ci.sp)

  SATE_95[irep, 1] = mean(colMeans(AFTbcf_fit1$tau)[true_label==1]) - mean(tau[true_label==1])
  SATE_95[irep, 2] = (mean(colMeans(AFTbcf_fit1$tau)[true_label==1]) - mean(tau[true_label==1]))^2
  l.ci.sp <- quantile(rowMeans(AFTbcf_fit1$tau[,true_label==1]),0.025)
  u.ci.sp <- quantile(rowMeans(AFTbcf_fit1$tau[,true_label==1]),0.975)
  SATE_95[irep, 3] = (mean(tau[true_label==1]) > l.ci.sp)*(mean(tau[true_label==1]) < u.ci.sp)
  SATE_95[irep, 4] = mean(u.ci.sp - l.ci.sp)

  SATE_99[irep, 1] = mean(colMeans(AFTbcf_fit2$tau)[true_label==1]) - mean(tau[true_label==1])
  SATE_99[irep, 2] = (mean(colMeans(AFTbcf_fit2$tau)[true_label==1]) - mean(tau[true_label==1]))^2
  l.ci.sp <- quantile(rowMeans(AFTbcf_fit2$tau[,true_label==1]),0.025)
  u.ci.sp <- quantile(rowMeans(AFTbcf_fit2$tau[,true_label==1]),0.975)
  SATE_99[irep, 3] = (mean(tau[true_label==1]) > l.ci.sp)*(mean(tau[true_label==1]) < u.ci.sp)
  SATE_99[irep, 4] = mean(u.ci.sp - l.ci.sp)
  
  SATE_obj[irep, 1] = mean(colMeans(obj$Theta)[true_label==1]) - mean(tau[true_label==1])
  SATE_obj[irep, 2] = (mean(colMeans(obj$Theta)[true_label==1]) - mean(tau[true_label==1]))^2
  l.ci.sp <- quantile(rowMeans(obj$Theta[,true_label==1]),0.025)
  u.ci.sp <- quantile(rowMeans(obj$Theta[,true_label==1]),0.975)
  SATE_obj[irep,3] = (mean(tau[true_label==1]) > l.ci.sp)*(mean(tau[true_label==1]) < u.ci.sp)
  SATE_obj[irep, 4] = mean(u.ci.sp - l.ci.sp)
  
  SATE_objsp[irep, 1] = mean(colMeans(obj_sp$Theta)[true_label==1]) - mean(tau[true_label==1])
  SATE_objsp[irep, 2] = (mean(colMeans(obj_sp$Theta)[true_label==1]) - mean(tau[true_label==1]))^2
  l.ci.sp <- quantile(rowMeans(obj_sp$Theta[,true_label==1]),0.025)
  u.ci.sp <- quantile(rowMeans(obj_sp$Theta[,true_label==1]),0.975)
  SATE_objsp[irep,3] = (mean(tau[true_label==1]) > l.ci.sp)*(mean(tau[true_label==1]) < u.ci.sp)
  SATE_objsp[irep, 4] = mean(u.ci.sp - l.ci.sp)

  s_1 = pnorm((AFTbcf_fit1$mu_trt - log(checktime))/AFTbcf_fit1$sigma)
  p_1 = (1-AFTbcf_fit1$pred_cureprob_1) + AFTbcf_fit1$pred_cureprob_1*s_1
  s_0 = pnorm((AFTbcf_fit1$mu_ctr - log(checktime))/AFTbcf_fit1$sigma)
  p_0 = (1-AFTbcf_fit1$pred_cureprob_0) + AFTbcf_fit1$pred_cureprob_0*s_0

  SurvPr_95[irep, 1] = mean(colMeans(p_1 - p_0) - true_survp)
  SurvPr_95[irep, 2] = sqrt(mean((colMeans(p_1 - p_0) - true_survp)^2))
  l.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.975))
  SurvPr_95[irep,3] = mean((true_survp > l.ci.sp) & (true_survp <= u.ci.sp))
  SurvPr_95[irep, 4] = mean(u.ci.sp - l.ci.sp)

  s_1 = pnorm((AFTbcf_fit2$mu_trt - log(checktime))/AFTbcf_fit2$sigma)
  p_1 = (1-AFTbcf_fit2$pred_cureprob_1) + AFTbcf_fit2$pred_cureprob_1*s_1
  s_0 = pnorm((AFTbcf_fit2$mu_ctr - log(checktime))/AFTbcf_fit2$sigma)
  p_0 = (1-AFTbcf_fit2$pred_cureprob_0) + AFTbcf_fit2$pred_cureprob_0*s_0

  SurvPr_99[irep, 1] = mean(colMeans(p_1 - p_0) - true_survp)
  SurvPr_99[irep, 2] = sqrt(mean((colMeans(p_1 - p_0) - true_survp)^2))
  l.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.975))
  SurvPr_99[irep,3] = mean((true_survp > l.ci.sp) & (true_survp <= u.ci.sp))
  SurvPr_99[irep, 4] = mean(u.ci.sp - l.ci.sp)
  
  mat0 = mat1 = obj$yhat.train
  mat1[, a==0] = mat1[,a==0] + obj$Theta[,a==0]
  p_1 = obj$mix.prop[,1]*pnorm((mat1 - log(checktime) - obj$locations[,1])/obj$sigma)
  mat0[, a==1] = mat0[,a==1] - obj$Theta[,a==1]
  p_0 = obj$mix.prop[,1]*pnorm((mat0 - log(checktime) - obj$locations[,1])/obj$sigma)
  for (i in 2:ncol(obj$mix.prop)) {
  p_1 = p_1 + obj$mix.prop[,i]*pnorm((mat1 - log(checktime) - obj$locations[,i])/obj$sigma)
  p_0 = p_0 + obj$mix.prop[,i]*pnorm((mat0 - log(checktime) - obj$locations[,i])/obj$sigma)
  }
  
  SurvPr_obj[irep, 1] = mean(colMeans(p_1 - p_0) - true_survp)
  SurvPr_obj[irep, 2] = sqrt(mean((colMeans(p_1 - p_0) - true_survp)^2))
  l.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.975))
  SurvPr_obj[irep,3] = mean((true_survp > l.ci.sp) & (true_survp <= u.ci.sp))
  SurvPr_obj[irep, 4] = mean(u.ci.sp - l.ci.sp)
  
  mat0 = mat1 = obj_sp$yhat.train
  mat1[, a==0] = mat1[,a==0] + obj_sp$Theta[,a==0]
  p_1 = pnorm((mat1 - log(checktime))/obj_sp$sigma)
  mat0[, a==1] = mat0[,a==1] - obj_sp$Theta[,a==1]
  p_0 = pnorm((mat0 - log(checktime))/obj_sp$sigma)

  SurvPr_objsp[irep, 1] = mean(colMeans(p_1 - p_0) - true_survp)
  SurvPr_objsp[irep, 2] = sqrt(mean((colMeans(p_1 - p_0) - true_survp)^2))
  l.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.025))
  u.ci.sp <- apply(p_1 - p_0, 2, function(x) quantile(x, prob=.975))
  SurvPr_objsp[irep,3] = mean((true_survp > l.ci.sp) & (true_survp <= u.ci.sp))
  SurvPr_objsp[irep, 4] = mean(u.ci.sp - l.ci.sp)
}

tmp = list('res_95' = res_95, 'res_99' = res_99, 'res_obj' = res_obj,'res_objsp' = res_objsp)
tmp2 = list('SATE_95' = SATE_95, 'SATE_99' = SATE_99, 'SATE_obj' = SATE_obj, 'SATE_objsp' = SATE_objsp)
tmp3 = list('SurvPr_95' = SurvPr_95, 'SurvPr_99' = SurvPr_99, 'SurvPr_obj' = SurvPr_obj, 'SurvPr_objsp' = SurvPr_objsp)
time.total <- proc.time() - time.begin
cat('Task complete in', time.total[3], 'seconds. Results saved in .csv.\n')
write.csv(tmp,paste(flag,'_res.csv',sep = ''))
write.csv(tmp2,paste(flag,'_SATE.csv',sep = ''))
write.csv(tmp3,paste(flag,'_SurvPr.csv',sep = ''))



