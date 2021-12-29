#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom RcppParallel RcppParallelLibs

Rcpp::loadModule(module = "TreeSamples", TRUE)

# Note: the code is established based on R package bcf @ https://github.com/jaredsmurray/bcf

.cp_quantile = function(x, num = 1000, cat_levels = 8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq == 1){
    ret = x[1]
    warning("A supplied covariate contains a single unique value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x), quantile(x, p=0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out = num)
    ret = q(ind)
  }
  return(ret)
}


.get_chain_tree_files = function(tree_path, chain_id){
  out <- list("con_trees" = paste0(tree_path, '/',"con_trees.", chain_id, ".txt"),
              "mod_trees" = paste0(tree_path, '/',"mod_trees.", chain_id, ".txt"),
              "cure_trees" = paste0(tree_path, '/',"cure_trees.", chain_id, ".txt"))
  return(out)
}

.get_do_type = function(n_cores){
  if(n_cores>1){
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    `%doType%`  <- foreach::`%dopar%`
  } else {
    cl <- NULL
    `%doType%`  <- foreach::`%do%`
  }

  do_type_config <- list('doType'  = `%doType%`,
                         'n_cores' = n_cores,
                         'cluster' = cl)

  return(do_type_config)
}

.cleanup_after_par = function(do_type_config){
  if(do_type_config$n_cores>1){
    parallel::stopCluster(do_type_config$cluster)
  }
}

.ident <- function(...){
  # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

# The following function is referenced from
# https://github.com/nchenderson/AFTrees
# density of inv-chisq
.dinvchisq <- function(x, nu){
  ans <- (nu/2)*log(nu/2) - lgamma(nu/2) - (nu/2 + 1)*log(x) - nu/(2*x)
  return(exp(ans))
}


.pchiconv <- function(p,nu){
  ## CDF for the random variable defined as nu/X + Z, X ~ chisq_nu, Z ~ N(1,1)
  ff <- function(v,x){
    ans <- pnorm(x-v, mean = 1, sd = 1)*.dinvchisq(v, nu)
    return(ans)
  }
  tau <- length(p)
  a <- rep(0, tau)
  for (k in 1:tau) {
    a[k] <- integrate(ff, x = p[k], lower = 0, upper = Inf)$value
  }
  return(a)
}

.qchiconv <- function(q, nu){
  gg <- function(x){
    .pchiconv(x, nu) - q
  }
  a <- uniroot(gg, interval = c(0,200))
  return(a$root)
}

.FindKappa <- function(q, sighat, nu){
  if(q < 0.905 & q > 0.895){
    Q <- 6.356677
  } else if (q < 0.995 & q > 0.985) {
    Q <- 27.17199
  } else if (q < .8 & q > .7){
    Q <- 3.858311
  } else if (q < .55 & q > .45){
    Q <- 2.531834
  } else if (q < .3 & q > .2){
    Q <- 1.569312
  } else {
    Q <- .qchiconv(q, nu)
  }
  Kappa <- sighat^2/Q
  return(Kappa)
}

#' @useDynLib TBAFTcure
#' @export
TBAFTcure <- function(y, status, a, x_control, x_moderate = x_control, x_cure,
                pihat, w = NULL, random_seed = sample.int(.Machine$integer.max,1),
                n_chains = 1, n_cores = n_chains,n_threads = max((RcppParallel::defaultNumThreads()-2)/n_cores, 1),
                nburn = 500, nsim = 1000, nthin = 1, update_interval = 100,
                ntree_control = 200,
                #sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                #sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3,
                ntree_cure = 200,
                sd_cure = 0.5,
                base_cure = 0.95,
                power_cure = 2,
                save_tree_directory = '.',
                nu = 3, lambda = NULL, sigq = 0.95, sighat = NULL,
                include_pi = "control", verbose = TRUE) {

  if(is.null(w)){
    w <- matrix(1,ncol = 1,nrow = length(y))
  }

  pihat = as.matrix(pihat)

  ### TODO range check on parameters

  x_c = matrix(x_control, ncol = ncol(x_control))
  x_m = matrix(x_moderate, ncol = ncol(x_moderate))
  x_cu = matrix(x_cure, ncol = ncol(x_cure))
  # include ps as a covariate
  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_c, pihat)
  }
  
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_m, pihat)
  }
  
  if(include_pi=="both" | include_pi == "cure" | include_pi == "control") {
    x_cu = cbind(x_cu, pihat)
  }

  ## set cutpoints
  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))
  cutpoint_list_cu = lapply(1:ncol(x_cu), function(i) .cp_quantile(x_cu[,i]))

  ## centering  log_y around mu_aft
  nullfit <- survreg(Surv(y, status) ~ 1, dist = 'lognormal')
  null_sig <- nullfit$scale
  null_inter <- nullfit$coefficients

  y_log_centered = log(y) - null_inter # logy - muhat_aft
  zeta = 4*null_sig # for variance hyper par of the leaf nodes
  #imr = dnorm(log(y), null_inter, null_sig)/pnorm(log(y), null_inter, null_sig, lower.tail = FALSE)
  #y_unobs = exp(status*log(y) + (1-status)*(null_sig + imr))
  #zeta_tmp = range(log(y_unobs))[2] - range(log(y_unobs))[1]

  if(is.null(sighat)){
    tempx <- cbind(x_c, x_m, a)
    tempx <- tempx[,!duplicated(t(tempx))]
    tempfit <- survreg(Surv(y,status) ~ tempx, dist = 'lognormal')
    sighat <- tempfit$scale # this is the sigma_hat_w in Henderson's work, for estiamting sigma
  }

  if(is.null(lambda)){
    qchi = qchisq(1.0 - sigq, nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  kappa <- .FindKappa(sigq, sighat, nu)

  dir = tempdir()
  perm = order(a, decreasing = TRUE) # so that the treated are displayed first

  cure_sd = 1.5
  
  RcppParallel::setThreadOptions(numThreads = n_threads)

  do_type_config <- .get_do_type(n_cores)
  `%doType%` <- do_type_config$doType


  chain_out <- foreach::foreach(iChain = 1:n_chains,
                                .export = c(".get_chain_tree_files","TBAFTcureRcpp")) %doType% {
    this_seed = random_seed + iChain - 1
    cat("Calling TBAFTcureRcpp From R \n")
    set.seed(this_seed)

    tree_files = .get_chain_tree_files(save_tree_directory, iChain)

    print(tree_files)

    fitbcf = TBAFTcureRcpp(y_ = y_log_centered[perm], a_ = a[perm], w_ = w[perm],
                                 x_con_ = t(x_c[perm,,drop=FALSE]), x_mod_ = t(x_m[perm,,drop=FALSE]), 
                                 x_cure_ = t(x_cu[perm,,drop = FALSE]), status_ = status[perm],
                                 x_con_info_list = cutpoint_list_c,
                                 x_mod_info_list = cutpoint_list_m,
                                 x_cure_info_list = cutpoint_list_cu,
                                 burn = nburn, nd = nsim, thin = nthin, ntree_cure = ntree_cure,
                                 ntree_mod = ntree_moderate, ntree_con = ntree_control,
                                 lambda = lambda, nu = nu,
                                 #con_sd = con_sd,
                                 #mod_sd = mod_sd,
                                 con_alpha = base_control,
                                 con_beta = power_control,
                                 mod_alpha = base_moderate,
                                 mod_beta = power_moderate,
                                 cure_alpha = base_cure,
                                 cure_beta = power_cure,
                                 cure_sd = cure_sd,
                                 kappa = kappa,
                                 sigest = sighat,
                                 zeta = zeta,
                                 treef_con_name_ = tree_files$con_trees,
                                 treef_mod_name_ = tree_files$mod_trees,
                                 treef_cure_name_ = tree_files$cure_trees,
                                 printevery = update_interval,
                                verbose_sigma = verbose)

    cat("TBAFTcureRcpp returned to R\n")
    
    gl <- fitbcf$glabel_post[, order(perm)] # stores the cured group labels
    curefit <- fitbcf$curefit_post[, order(perm)]

    ac = fitbcf$m_post[, order(perm)] # stores allfit_con

    Tm = fitbcf$b_post[, order(perm)]*(1.0/(fitbcf$b1 - fitbcf$b0)) # stores the allfit_mod/bscale, the pure tree fit

    Tc = ac*(1.0/fitbcf$msd) #stores allfit_con/mscale, the pure tree fit

    tau_post = fitbcf$b_post[, order(perm)]

    mu_post = null_inter + Tc*fitbcf$msd + Tm*fitbcf$b0 #stores E[Y|0,X] for everyone, note: without intercept

    mu_post_trt = null_inter + Tc*fitbcf$msd + Tm*fitbcf$b1 # stores E[Y|1,X] for everyone note: without intercept

    varcnt_con_post = fitbcf$varcnt_con_post
    varcnt_mod_post = fitbcf$varcnt_mod_post
    varcnt_cure_post = fitbcf$varcnt_cure_post
    pred_cureprob_0_post = fitbcf$pred_cureprob_0_post
    pred_cureprob_1_post = fitbcf$pred_cureprob_1_post


    #chain_out = list()
    #chain_out[[1]] <-
    list(sigma = fitbcf$sigma,
         yhat = null_inter + fitbcf$yhat_post[,order(perm)],
         muy = null_inter,
         mu_ctr  = mu_post,
         mu_trt = mu_post_trt,
         tau = tau_post,
         mu_scale = fitbcf$msd,
         tau_scale = fitbcf$bsd,
         b0 = fitbcf$b0,
         b1 = fitbcf$b1,
         perm = perm,
         include_pi = include_pi,
         random_seed=this_seed,
         varcnt_con = varcnt_con_post,
         varcnt_mod = varcnt_mod_post,
         varcnt_cure = varcnt_cure_post,
         cured_label = gl,
         cured_fit = curefit,
         pred_cureprob_0 = pred_cureprob_0_post,
         pred_cureprob_1 = pred_cureprob_1_post
    )

  }

  all_sigma = c()
  all_mu_scale = c()
  all_tau_scale = c()

  all_b0 = c()
  all_b1 = c()

  all_yhat = c()
  all_mu_ctr   = c()
  all_mu_trt   = c()
  all_tau  = c()

  all_varcnt_con = c()
  all_varcnt_mod = c()
  all_varcnt_cure = c()

  all_cured_label = c()
  all_cured_fit = c()
  all_pred_cureprob_0 = c()
  all_pred_cureprob_1 = c()


  chain_list=list()

  #n_iter = length(chain_out[[1]]$sigma)
  #
  for (iChain in 1:n_chains){
    sigma <- chain_out[[iChain]]$sigma
    mu_scale  <- chain_out[[iChain]]$mu_scale
    tau_scale <- chain_out[[iChain]]$tau_scale

    b0 <- chain_out[[iChain]]$b0
    b1 <- chain_out[[iChain]]$b1

    yhat  <- chain_out[[iChain]]$yhat
    tau <- chain_out[[iChain]]$tau
    mu_ctr  <- chain_out[[iChain]]$mu_ctr
    mu_trt  <- chain_out[[iChain]]$mu_trt

    varcnt_con <- chain_out[[iChain]]$varcnt_con
    varcnt_mod <- chain_out[[iChain]]$varcnt_mod
    varcnt_cure <- chain_out[[iChain]]$varcnt_cure
    
    cured_label <- chain_out[[iChain]]$cured_label
    cured_fit <- chain_out[[iChain]]$cured_fit
    pred_cureprob_0 <- chain_out[[iChain]]$pred_cureprob_0
    pred_cureprob_1 <- chain_out[[iChain]]$pred_cureprob_1

    # -----------------------------
    # Support Old Output
    # -----------------------------
    all_sigma = c(all_sigma, sigma)
    all_mu_scale = c(all_mu_scale,  mu_scale)
    all_tau_scale = c(all_tau_scale, tau_scale)
    all_b0 = c(all_b0, b0)
    all_b1 = c(all_b1, b1)

    all_yhat = rbind(all_yhat, yhat)
    all_mu_ctr  = rbind(all_mu_ctr, mu_ctr)
    all_mu_trt  = rbind(all_mu_trt, mu_trt)
    all_tau = rbind(all_tau, tau)

    all_varcnt_con = rbind(all_varcnt_con,varcnt_con)
    all_varcnt_mod = rbind(all_varcnt_mod,varcnt_mod)
    all_varcnt_cure = rbind(all_varcnt_cure, varcnt_cure)
    
    all_cured_label = rbind(all_cured_label, cured_label)
    all_cured_fit = rbind(all_cured_fit, cured_fit)
    all_pred_cureprob_0 = rbind(all_pred_cureprob_0, pred_cureprob_0)
    all_pred_cureprob_1 = rbind(all_pred_cureprob_1, pred_cureprob_1)
    


    # -----------------------------
    # Make the MCMC Object
    # -----------------------------

    scalar_df <- data.frame("sigma" = sigma,
                            "tau_bar" = matrixStats::rowWeightedMeans(tau, w),
                            "mu_ctr_bar"  = matrixStats::rowWeightedMeans(mu_ctr, w),
                            "mu_trt_bar"  = matrixStats::rowWeightedMeans(mu_trt, w),
                            "yhat_bar" = matrixStats::rowWeightedMeans(yhat, w),
                            "mu_scale" = mu_scale,
                            "tau_scale" = tau_scale,
                            "b0"  = b0,
                            "b1"  = b1,
                            "cured_fit_bar" = matrixStats::rowWeightedMeans(pnorm(cured_fit), w),
                            "cured_labels_bar" = cured_label,
                            "pred_cureprob_0_bar"  = matrixStats::rowWeightedMeans(pred_cureprob_0, w),
                            "pred_cureprob_1_bar"  = matrixStats::rowWeightedMeans(pred_cureprob_1, w))
    



    chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    # -----------------------------
    # Sanity Check Constants Accross Chains
    # -----------------------------
    #if(chain_out[[iChain]]$con_sd   != chain_out[[1]]$con_sd)     stop("con_sd not consistent between chains for no reason")
    #if(chain_out[[iChain]]$mod_sd   != chain_out[[1]]$mod_sd)     stop("mod_sd not consistent between chains for no reason")
    #if(chain_out[[iChain]]$muy      != chain_out[[1]]$muy)        stop("muy not consistent between chains for no reason")
    if(chain_out[[iChain]]$include_pi != chain_out[[1]]$include_pi) stop("include_pi not consistent between chains for no reason")
    if(any(chain_out[[iChain]]$perm   != chain_out[[1]]$perm))      stop("perm not consistent between chains for no reason")
  }

  fitObj <- list(sigma = all_sigma,
                 yhat = all_yhat,
                 #sdy = chain_out[[1]]$sdy,
                 muy = chain_out[[1]]$muy,
                 mu_ctr  = all_mu_ctr,
                 mu_trt  = all_mu_trt,
                 tau = all_tau,
                 mu_scale = all_mu_scale,
                 tau_scale = all_tau_scale,
                 b0 = all_b0,
                 b1 = all_b1,
                 perm = perm,
                 include_pi = chain_out[[1]]$include_pi,
                 random_seed = chain_out[[1]]$random_seed,
                 coda_chains = coda::as.mcmc.list(chain_list),
                 raw_chains = chain_out,
                 varcnt_con = all_varcnt_con,
                 varcnt_mod = all_varcnt_mod,
                 varcnt_cure = all_varcnt_cure,
                 cured_label = all_cured_label,
                 cured_fit = all_cured_fit,
                 pred_cureprob_0 = all_pred_cureprob_0,
                 pred_cureprob_1 = all_pred_cureprob_1,
                 coda_chains = coda::as.mcmc.list(chain_list),
                 raw_chains = chain_out)

  attr(fitObj, "class") <- "tbaftcure"

  .cleanup_after_par(do_type_config)

  return(fitObj)
}

#' @export predict.tbaftcure
#' @export
predict.tbaftcure <- function(object, x_pred_control, x_pred_moderate, x_pred_cure, 
                           pihat_pred, a_pred,
                           save_tree_directory, ncores = 1,
                           ...) {
  if(any(is.na(x_pred_control))) stop("Missing values in x_pred_control")
  if(any(is.na(x_pred_moderate))) stop("Missing values in x_pred_moderate")
  if(any(is.na(x_pred_cure))) stop("Missing values in x_pred_cure")
  if(any(is.na(pihat_pred))) stop("Missing values in pihat_pred")
  if(any(is.na(a_pred))) stop("Missing values in a_pred")

  pihat_pred = as.matrix(pihat_pred)

  x_pc = matrix(x_pred_control, ncol = ncol(x_pred_control))
  x_pm = matrix(x_pred_moderate, ncol = ncol(x_pred_moderate))
  x_pcu = matrix(x_pred_cure, ncol = ncol(x_pred_cure))

  if(object$include_pi == "both" | object$include_pi == "control") {
    x_pc = cbind(x_pc, pihat_pred)
  }

  if(object$include_pi == "both" | object$include_pi == "moderate") {
    x_pm = cbind(x_pm, pihat_pred)
  }
  
  if(object$include_pi == "both" | object$include_pi == "control" | object$include_pi == "cure") {
    x_pcu = cbind(x_pcu, pihat_pred)
  }

  cat("Starting Prediction \n")
  n_chains = length(object$coda_chains)

  do_type_config <- .get_do_type(ncores)
  `%doType%` <- do_type_config$doType

  chain_out <- foreach::foreach(iChain = 1:n_chains) %doType% {
    tree_files = .get_chain_tree_files(save_tree_directory, iChain)
    cat("Starting to Predict Chain ", iChain, "\n")
    # Note: what obtained below are prediction made from the saved ndraw posterior draws of trees
    modts = TreeSamples$new()
    modts$load(tree_files$mod_trees)
    Tm = modts$predict(t(x_pm)) #pure mod tree fit, dim(ndraws,n)

    conts = TreeSamples$new()
    conts$load(tree_files$con_trees)
    Tc = conts$predict(t(x_pc)) # pure con tree fit
    
    curets = TreeSamples$new()
    curets$load(tree_files$cure_trees)
    Tcure = curets$predict(t(x_pcu))

    list(Tm = Tm, Tc = Tc, Tcure = Tcure)
  }

  all_yhat = c()
  all_mu = c()
  all_tau = c()
  
  all_cure = c()

  chain_list = list()

  for (iChain in 1:n_chains) {
    Tm = chain_out[[iChain]]$Tm
    Tc = chain_out[[iChain]]$Tc
    Tcure = chain_out[[iChain]]$Tcure

    this_chain_bcf_out = object$raw_chains[[iChain]]

    null_inter = object$muy
    b1 = this_chain_bcf_out$b1
    b0 = this_chain_bcf_out$b0
    mu_scale = this_chain_bcf_out$mu_scale

    # get tau, mu, and y: the three parts that we care about in prediction
    mu = null_inter + Tc*mu_scale + Tm*b0
    tau = (b1 - b0)*Tm #mu+tau = mu_trt
    yhat = mu + t(t(tau)*a_pred)

    all_yhat = rbind(all_yhat, yhat)
    all_mu = rbind(all_mu, mu)
    all_tau = rbind(all_tau, tau)
    
    #cured_intc = object$cured_intc
    cured_fit = Tcure# + cured_intc ## BUT mind the dim, can we directly add them up (matrix+ vec)?
    all_cure = rbind(all_cure, pnorm(cured_fit)) # the posterior probability of being uncured
    
    scalar_df <- data.frame("tau_bar" = matrixStats::rowWeightedMeans(tau, w = NULL),
                            "mu_bar" = matrixStats::rowWeightedMeans(mu, w = NULL),
                            "yhat_bar" = matrixStats::rowWeightedMeans(yhat, w = NULL),
                            "uncured_prob_bar" = matrixStats::rowWeightedMeans(pnorm(cured_fit), w = NULL))

    chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
  }

  .cleanup_after_par(do_type_config)

  list(tau = all_tau,
       mu = all_mu,
       yhat = all_yhat,
       uncured_prob = all_cure,
       coda_chains = coda::as.mcmc.list(chain_list))
}
