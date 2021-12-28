// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// TBAFTcureRcpp
List TBAFTcureRcpp(NumericVector y_, NumericVector a_, NumericVector w_, NumericVector x_con_, NumericVector x_mod_, NumericVector x_cure_, IntegerVector status_, List x_con_info_list, List x_mod_info_list, List x_cure_info_list, int burn, int nd, int thin, int ntree_cure, int ntree_mod, int ntree_con, double lambda, double nu, double con_alpha, double con_beta, double mod_alpha, double mod_beta, double cure_alpha, double cure_beta, double cure_sd, double kappa, double sigest, double zeta, CharacterVector treef_con_name_, CharacterVector treef_mod_name_, CharacterVector treef_cure_name_, int printevery, double trt_init, bool verbose_sigma);
static SEXP _TBAFTcure_TBAFTcureRcpp_try(SEXP y_SEXP, SEXP a_SEXP, SEXP w_SEXP, SEXP x_con_SEXP, SEXP x_mod_SEXP, SEXP x_cure_SEXP, SEXP status_SEXP, SEXP x_con_info_listSEXP, SEXP x_mod_info_listSEXP, SEXP x_cure_info_listSEXP, SEXP burnSEXP, SEXP ndSEXP, SEXP thinSEXP, SEXP ntree_cureSEXP, SEXP ntree_modSEXP, SEXP ntree_conSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP con_alphaSEXP, SEXP con_betaSEXP, SEXP mod_alphaSEXP, SEXP mod_betaSEXP, SEXP cure_alphaSEXP, SEXP cure_betaSEXP, SEXP cure_sdSEXP, SEXP kappaSEXP, SEXP sigestSEXP, SEXP zetaSEXP, SEXP treef_con_name_SEXP, SEXP treef_mod_name_SEXP, SEXP treef_cure_name_SEXP, SEXP printeverySEXP, SEXP trt_initSEXP, SEXP verbose_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a_(a_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_con_(x_con_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_mod_(x_mod_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_cure_(x_cure_SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type status_(status_SEXP);
    Rcpp::traits::input_parameter< List >::type x_con_info_list(x_con_info_listSEXP);
    Rcpp::traits::input_parameter< List >::type x_mod_info_list(x_mod_info_listSEXP);
    Rcpp::traits::input_parameter< List >::type x_cure_info_list(x_cure_info_listSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type ntree_cure(ntree_cureSEXP);
    Rcpp::traits::input_parameter< int >::type ntree_mod(ntree_modSEXP);
    Rcpp::traits::input_parameter< int >::type ntree_con(ntree_conSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type con_alpha(con_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type con_beta(con_betaSEXP);
    Rcpp::traits::input_parameter< double >::type mod_alpha(mod_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mod_beta(mod_betaSEXP);
    Rcpp::traits::input_parameter< double >::type cure_alpha(cure_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type cure_beta(cure_betaSEXP);
    Rcpp::traits::input_parameter< double >::type cure_sd(cure_sdSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type sigest(sigestSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type treef_con_name_(treef_con_name_SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type treef_mod_name_(treef_mod_name_SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type treef_cure_name_(treef_cure_name_SEXP);
    Rcpp::traits::input_parameter< int >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< double >::type trt_init(trt_initSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_sigma(verbose_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(TBAFTcureRcpp(y_, a_, w_, x_con_, x_mod_, x_cure_, status_, x_con_info_list, x_mod_info_list, x_cure_info_list, burn, nd, thin, ntree_cure, ntree_mod, ntree_con, lambda, nu, con_alpha, con_beta, mod_alpha, mod_beta, cure_alpha, cure_beta, cure_sd, kappa, sigest, zeta, treef_con_name_, treef_mod_name_, treef_cure_name_, printevery, trt_init, verbose_sigma));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _TBAFTcure_TBAFTcureRcpp(SEXP y_SEXP, SEXP a_SEXP, SEXP w_SEXP, SEXP x_con_SEXP, SEXP x_mod_SEXP, SEXP x_cure_SEXP, SEXP status_SEXP, SEXP x_con_info_listSEXP, SEXP x_mod_info_listSEXP, SEXP x_cure_info_listSEXP, SEXP burnSEXP, SEXP ndSEXP, SEXP thinSEXP, SEXP ntree_cureSEXP, SEXP ntree_modSEXP, SEXP ntree_conSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP con_alphaSEXP, SEXP con_betaSEXP, SEXP mod_alphaSEXP, SEXP mod_betaSEXP, SEXP cure_alphaSEXP, SEXP cure_betaSEXP, SEXP cure_sdSEXP, SEXP kappaSEXP, SEXP sigestSEXP, SEXP zetaSEXP, SEXP treef_con_name_SEXP, SEXP treef_mod_name_SEXP, SEXP treef_cure_name_SEXP, SEXP printeverySEXP, SEXP trt_initSEXP, SEXP verbose_sigmaSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_TBAFTcure_TBAFTcureRcpp_try(y_SEXP, a_SEXP, w_SEXP, x_con_SEXP, x_mod_SEXP, x_cure_SEXP, status_SEXP, x_con_info_listSEXP, x_mod_info_listSEXP, x_cure_info_listSEXP, burnSEXP, ndSEXP, thinSEXP, ntree_cureSEXP, ntree_modSEXP, ntree_conSEXP, lambdaSEXP, nuSEXP, con_alphaSEXP, con_betaSEXP, mod_alphaSEXP, mod_betaSEXP, cure_alphaSEXP, cure_betaSEXP, cure_sdSEXP, kappaSEXP, sigestSEXP, zetaSEXP, treef_con_name_SEXP, treef_mod_name_SEXP, treef_cure_name_SEXP, printeverySEXP, trt_initSEXP, verbose_sigmaSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _TBAFTcure_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("List(*TBAFTcureRcpp)(NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,IntegerVector,List,List,List,int,int,int,int,int,int,double,double,double,double,double,double,double,double,double,double,double,double,CharacterVector,CharacterVector,CharacterVector,int,double,bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _TBAFTcure_RcppExport_registerCCallable() { 
    R_RegisterCCallable("TBAFTcure", "_TBAFTcure_TBAFTcureRcpp", (DL_FUNC)_TBAFTcure_TBAFTcureRcpp_try);
    R_RegisterCCallable("TBAFTcure", "_TBAFTcure_RcppExport_validate", (DL_FUNC)_TBAFTcure_RcppExport_validate);
    return R_NilValue;
}

RcppExport SEXP _rcpp_module_boot_TreeSamples();

static const R_CallMethodDef CallEntries[] = {
    {"_TBAFTcure_TBAFTcureRcpp", (DL_FUNC) &_TBAFTcure_TBAFTcureRcpp, 34},
    {"_rcpp_module_boot_TreeSamples", (DL_FUNC) &_rcpp_module_boot_TreeSamples, 0},
    {"_TBAFTcure_RcppExport_registerCCallable", (DL_FUNC) &_TBAFTcure_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_TBAFTcure(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
