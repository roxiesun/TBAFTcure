# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

TBAFTcureRcpp <- function(y_, a_, w_, x_con_, x_mod_, x_cure_, status_, x_con_info_list, x_mod_info_list, x_cure_info_list, burn, nd, thin, ntree_cure, ntree_mod, ntree_con, lambda, nu, con_alpha, con_beta, mod_alpha, mod_beta, cure_alpha, cure_beta, cure_sd, kappa, sigest, zeta, treef_con_name_, treef_mod_name_, treef_cure_name_, printevery = 100L, trt_init = 1.0, verbose_sigma = FALSE) {
    .Call('_TBAFTcure_TBAFTcureRcpp', PACKAGE = 'TBAFTcure', y_, a_, w_, x_con_, x_mod_, x_cure_, status_, x_con_info_list, x_mod_info_list, x_cure_info_list, burn, nd, thin, ntree_cure, ntree_mod, ntree_con, lambda, nu, con_alpha, con_beta, mod_alpha, mod_beta, cure_alpha, cure_beta, cure_sd, kappa, sigest, zeta, treef_con_name_, treef_mod_name_, treef_cure_name_, printevery, trt_init, verbose_sigma)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_TBAFTcure_RcppExport_registerCCallable', PACKAGE = 'TBAFTcure')
})
