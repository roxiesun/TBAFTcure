// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_TBAFTcure_RCPPEXPORTS_H_GEN_
#define RCPP_TBAFTcure_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace TBAFTcure {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("TBAFTcure", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("TBAFTcure", "_TBAFTcure_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in TBAFTcure");
            }
        }
    }

    inline List TBAFTcureRcpp(NumericVector y_, NumericVector a_, NumericVector w_, NumericVector x_con_, NumericVector x_mod_, NumericVector x_cure_, IntegerVector status_, List x_con_info_list, List x_mod_info_list, List x_cure_info_list, int burn, int nd, int thin, int ntree_cure, int ntree_mod, int ntree_con, double lambda, double nu, double con_alpha, double con_beta, double mod_alpha, double mod_beta, double cure_alpha, double cure_beta, double cure_sd, double kappa, double sigest, double zeta, CharacterVector treef_con_name_, CharacterVector treef_mod_name_, CharacterVector treef_cure_name_, int printevery = 100, double trt_init = 1.0, bool verbose_sigma = false) {
        typedef SEXP(*Ptr_TBAFTcureRcpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_TBAFTcureRcpp p_TBAFTcureRcpp = NULL;
        if (p_TBAFTcureRcpp == NULL) {
            validateSignature("List(*TBAFTcureRcpp)(NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,NumericVector,IntegerVector,List,List,List,int,int,int,int,int,int,double,double,double,double,double,double,double,double,double,double,double,double,CharacterVector,CharacterVector,CharacterVector,int,double,bool)");
            p_TBAFTcureRcpp = (Ptr_TBAFTcureRcpp)R_GetCCallable("TBAFTcure", "_TBAFTcure_TBAFTcureRcpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_TBAFTcureRcpp(Shield<SEXP>(Rcpp::wrap(y_)), Shield<SEXP>(Rcpp::wrap(a_)), Shield<SEXP>(Rcpp::wrap(w_)), Shield<SEXP>(Rcpp::wrap(x_con_)), Shield<SEXP>(Rcpp::wrap(x_mod_)), Shield<SEXP>(Rcpp::wrap(x_cure_)), Shield<SEXP>(Rcpp::wrap(status_)), Shield<SEXP>(Rcpp::wrap(x_con_info_list)), Shield<SEXP>(Rcpp::wrap(x_mod_info_list)), Shield<SEXP>(Rcpp::wrap(x_cure_info_list)), Shield<SEXP>(Rcpp::wrap(burn)), Shield<SEXP>(Rcpp::wrap(nd)), Shield<SEXP>(Rcpp::wrap(thin)), Shield<SEXP>(Rcpp::wrap(ntree_cure)), Shield<SEXP>(Rcpp::wrap(ntree_mod)), Shield<SEXP>(Rcpp::wrap(ntree_con)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(nu)), Shield<SEXP>(Rcpp::wrap(con_alpha)), Shield<SEXP>(Rcpp::wrap(con_beta)), Shield<SEXP>(Rcpp::wrap(mod_alpha)), Shield<SEXP>(Rcpp::wrap(mod_beta)), Shield<SEXP>(Rcpp::wrap(cure_alpha)), Shield<SEXP>(Rcpp::wrap(cure_beta)), Shield<SEXP>(Rcpp::wrap(cure_sd)), Shield<SEXP>(Rcpp::wrap(kappa)), Shield<SEXP>(Rcpp::wrap(sigest)), Shield<SEXP>(Rcpp::wrap(zeta)), Shield<SEXP>(Rcpp::wrap(treef_con_name_)), Shield<SEXP>(Rcpp::wrap(treef_mod_name_)), Shield<SEXP>(Rcpp::wrap(treef_cure_name_)), Shield<SEXP>(Rcpp::wrap(printevery)), Shield<SEXP>(Rcpp::wrap(trt_init)), Shield<SEXP>(Rcpp::wrap(verbose_sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_TBAFTcure_RCPPEXPORTS_H_GEN_