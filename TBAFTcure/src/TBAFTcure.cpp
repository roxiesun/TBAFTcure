#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rand_gen.h"
#include "tree.h"
#include "info.h"
#include "funcs.h"
#include "bd.h"
#include "logging.h"
#include "changerule.h"
#include "swaprule.h"

using namespace Rcpp;


/* model:
 cured fraction: probit Pr(G_i = 1) = \Phi(g_c(x,A)
 for uncured G_i = 1, log(t) = y = mu(x) + b(x)*a  + w
 w ~ N(0, \sigma^2)
 status = 1, observed; 0, censored;
 If censored, ylat ~ TN(mu(x) + b(x)a, \sigma^2, yobs). And ylat should be used in updaing the trees and labels.
 */


//[[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]

List TBAFTcureRcpp(NumericVector y_, NumericVector a_, NumericVector w_, NumericVector x_con_, NumericVector x_mod_, NumericVector x_cure_,  IntegerVector status_, List x_con_info_list, List x_mod_info_list, List x_cure_info_list, int burn, int nd, int thin, int ntree_cure, int ntree_mod, int ntree_con, double lambda, double nu, double con_alpha, double con_beta, double mod_alpha, double mod_beta, double cure_alpha, double cure_beta, double cure_sd, double kappa, double sigest, double zeta, CharacterVector treef_con_name_, CharacterVector treef_mod_name_, CharacterVector treef_cure_name_, int printevery=100,  double trt_init=1.0, bool verbose_sigma=false)
{

    std::string treef_con_name = as<std::string>(treef_con_name_);
    std::ofstream treef_con(treef_con_name.c_str());

    std::string treef_mod_name = as<std::string>(treef_mod_name_);
    std::ofstream treef_mod(treef_mod_name.c_str());
    
    std::string treef_cure_name = as<std::string>(treef_cure_name_);
    std::ofstream treef_cure(treef_cure_name.c_str());

    
    RNGScope scope;
    RNG gen; // one random number generator used in all draws
    
    
    Logger logger = Logger();
    char logBuff[100];

    bool log_level = false;

    logger.setLevel(log_level);
    logger.log("=================================");
    logger.log("Starting up TBAFTcure: ");
    logger.log("=================================");
    if (log_level) {
        logger.getVectorHead(y_, logBuff);
        Rcout << "y: " << logBuff << "\n";
        logger.getVectorHead(a_, logBuff);
        Rcout << "a: " << logBuff << "\n";
        logger.getVectorHead(w_, logBuff);
        Rcout << "w: " << logBuff << "\n";
    }

    //logger.log("Algorithm is Weighted");
    logger.log("");

    /* **** Read format y **** */
    std::vector<double> yobs, y;
    std::vector<int> delta, glabel;// y is for true survival time, yobs is the observed event/censoring time, delta is the event indicator
    double miny = INFINITY, maxy = -INFINITY;
    sinfo allys; // suff stat for all y, used to initialized the barts
    int num_censored = 0;
    
    for (NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
        yobs.push_back(*it);
        y.push_back(*it);
        if (*it<miny) miny=*it;
        if (*it>maxy) maxy=*it;
        allys.sy += *it; //\sum y
    }
    
    // read in status: 1 for event, 0 for censored
    for (IntegerVector::iterator it = status_.begin(); it!=status_.end(); ++it) {
        delta.push_back(*it);
        num_censored += (1 - *it);
        glabel.push_back(1);
    }
    
    Rcout << glabel[1] << "," << glabel[3] << endl;

    size_t n = y.size();
    allys.n = n;
    /*int* glabel = new int[n];
    for (size_t k=0; k<n; k++) glabel[k] = delta[k]; // initialized glabel*/
    double ybar = allys.sy/n; // sample mean
    // introduce the latent variable for cured group label
    sinfo allgs;
    std::vector<double> g;
    double lat = 0.;
    for (size_t k=0; k<n; k++) {
        if (glabel[k]>0) {
            lat = rtnormlo0(0);
            g.push_back(lat);
        } else {
            lat = -rtnormlo0(0);
            g.push_back(lat);
        }
        allgs.sy += lat;
    }
    allgs.n = g.size();
    if (allgs.n != n) {
        stop("dimension mismatch for the cured latent.");
    }
    
    

    /* **** Read format weights **** */
    double* w = new double[n];

    for (int j=0; j<n; j++) {
        w[j] = w_[j];
    }

    /* **** Read format x_con **** */
    std::vector<double> x_con;
    for (NumericVector::iterator it=x_con_.begin(); it!=x_con_.end(); ++it) {
        x_con.push_back(*it);
    }
    size_t p_con = x_con.size()/n;
    std::vector<size_t> nv_con(p_con,0.);

    Rcout << "Using " << p_con << " control variables." << std::endl;

    // cutpoints;
    xinfo xi_con;
    xi_con.resize(p_con);
    for (int i=0; i<p_con; ++i) {
        NumericVector tmp = x_con_info_list[i];
        std::vector<double> tmp2;
        for (size_t j=0; j<tmp.size(); ++j) {
            tmp2.push_back(tmp[j]);
        }
        xi_con[i] = tmp2;
        tmp2.clear();
    }

    /* **** Read format x_mod **** */
    int ntrt = 0;
    for (size_t i=0; i<n; ++i) {
        if (a_[i]>0) ntrt +=1; // Note: binary treatment!!
    }
    std::vector<double> x_mod;
    for (NumericVector::iterator it=x_mod_.begin(); it!=x_mod_.end(); ++it) {
        x_mod.push_back(*it);
    }
    size_t p_mod = x_mod.size()/n;
    std::vector<size_t> nv_mod(p_mod,0.);

    Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;

    // cutpoints
    xinfo xi_mod;

    xi_mod.resize(p_mod);
    for (int i=0; i<p_mod; ++i) {
        NumericVector tmp = x_mod_info_list[i];
        std::vector<double> tmp2;
        for (size_t j=0; j<tmp.size(); ++j) {
            tmp2.push_back(tmp[j]);
        }
        xi_mod[i] = tmp2;
        tmp2.clear();
    }
    
    /* **** Read format x_cure **** */
    std::vector<double> x_cure;
    for (NumericVector::iterator it=x_cure_.begin(); it!=x_cure_.end(); ++it) {
        x_cure.push_back(*it);
    }
    size_t p_cure = x_cure.size()/n;
    std::vector<size_t> nv_cure(p_cure,0.);

    Rcout << "Using " << p_cure << " predictors for the cure fraction." << std::endl;

    // cutpoints;
    xinfo xi_cure;
    xi_cure.resize(p_cure);
    for (int i=0; i<p_cure; ++i) {
        NumericVector tmp = x_cure_info_list[i];
        std::vector<double> tmp2;
        for (size_t j=0; j<tmp.size(); ++j) {
            tmp2.push_back(tmp[j]);
        }
        xi_cure[i] = tmp2;
        tmp2.clear();
    }
    
    // Read format predxcure_1, predxcure_0
    std::vector<double> pred_x_cure_1;
    std::vector<double> pred_x_cure_0;
    for (NumericVector::iterator it= x_cure_.begin(); it!=x_cure_.end(); ++it) {
        pred_x_cure_0.push_back(*it);
        pred_x_cure_1.push_back(*it);
    }
    for (size_t k=0; k<n; k++) {
        pred_x_cure_1[(k+1)*p_cure - 2] = 1.;
        pred_x_cure_0[(k+1)*p_cure - 2] = 0.;
    }
    
    
    /* **** Set up the model **** */

    // trees
    std::vector<tree> t_mod(ntree_mod);
    for (size_t i=0; i<ntree_mod; i++) {
        t_mod[i].setmu(trt_init/(double)ntree_mod);
    }

    std::vector<tree> t_con(ntree_con);
    for (size_t i=0; i<ntree_con; i++) {
        t_con[i].setmu(ybar/(double)ntree_con);
    }
    
    std::vector<tree> t_cure(ntree_cure);
    for (size_t i=0; i<ntree_cure; i++) {
        t_cure[i].setmu(allgs.sy/((double)(ntree_cure*n)));
    }

    // prior
    // scale parameters for b, the modifier part
    //double bscale_prec = 2;
    double bscale0 = -1.0; // are the two bscales ps ?
    double bscale1 = 1.0;

    //double mscale_prec = 1.0; // half cauchy prior for muscale
    double mscale = 1.0;
    //double delta_con = 1.0; // half normal prior for tau scale, the homo treatment effect
    //double delta_mod = 1.0;

    pinfo pi_mod;
    pi_mod.pbd = 0.5;
    pi_mod.pb = 0.25;
    pi_mod.pchange = 0.4;
    pi_mod.pswap = 0.1;
    

    pi_mod.alpha = mod_alpha;
    pi_mod.beta = mod_beta;
    pi_mod.tau = zeta/(8*sqrt(ntree_mod));  // sigma_mu, variance on leaf parameters, =zeta/(2k*sqrt(ntree))
    pi_mod.sigma = sigest; // residual variance is \sigma^2_y / bscale^2

    pinfo pi_con;
    pi_con.pbd = 0.5;
    pi_con.pb = 0.25;
    pi_con.pchange = 0.4;
    pi_con.pswap = 0.1;

    pi_con.alpha = con_alpha;
    pi_con.beta = con_beta;
    pi_con.tau = zeta/(4*sqrt(ntree_con)); //con_sd/(sqrt(delta_con)*sqrt((double) ntree_con));
    pi_con.sigma = sigest/fabs(mscale);

    double sigma = sigest;
    
    pinfo pi_cure;
    pi_cure.pbd = 0.5;
    pi_cure.pb = 0.25;
    pi_cure.pchange = 0.4;
    pi_cure.pswap = 0.1;
    
    pi_cure.alpha = cure_alpha;
    pi_cure.beta = cure_beta;
    pi_cure.tau = cure_sd/sqrt((double) ntree_cure);
    pi_cure.sigma = 1.0;
    


    // Initialize dinfo
    double* allfit_con = new double[n]; // sum of fit of all trees
    for (size_t i=0; i<n; i++) allfit_con[i] = ybar;
    double* r_con = new double[n]; // residual for tree_con
    dinfo di_con;
    di_con.N = n;
    di_con.p = p_con;
    di_con.x = &x_con[0];
    di_con.y = r_con; //Note: y for each draw is the residual!! y - indiv_location - (allfit - ftemp[fit of current tree]) = y - tau_Li - allfit + ftemp

    // dinfo for treatment effect funtion b(x) or tau(x);
    double* allfit_mod = new double[n];
    for (size_t i=0; i<n; i++) allfit_mod[i] = (a_[i]*bscale1 + (1-a_[i])*bscale0)*trt_init;
    double* r_mod = new double[n]; // residual
    dinfo di_mod;
    di_mod.N = n;
    di_mod.p = p_mod;
    di_mod.x = &x_mod[0];
    di_mod.y = r_mod;
                       
    
    //dinfo for prediction cure data
    dinfo di_pred_cure_1;
    di_pred_cure_1.N = n;
    di_pred_cure_1.p = p_cure;
    di_pred_cure_1.x = &pred_x_cure_1[0];
    
    dinfo di_pred_cure_0;
    di_pred_cure_0.N = n;
    di_pred_cure_0.p = p_cure;
    di_pred_cure_0.x = &pred_x_cure_0[0];
    
    
    //fitpred_cure_1, fitpred_cure_0 stores tree-fitted values for the test data
    double* fitpred_cure_0 = new double[n];
    double* fitpred_cure_1 = new double[n];
                       
    /*dinfo for the cured model */
    double* allfit_cure = new double[n]; double* cure_intc = new double [1]; cure_intc[0] = 0.;
    for (size_t i=0; i<n; i++) allfit_cure[i] = allgs.sy/n;
    double* r_cure = new double[n];
    dinfo di_cure;
    di_cure.N = n;
    di_cure.p = p_cure;
    di_cure.x = &x_cure[0];
    di_cure.y = r_cure;
                       
    
    // ------------------------------
    // store the fits
    double* allfit = new double[n]; //yhat
    for (size_t i=0; i<n; i++) {
        allfit[i] = allfit_mod[i] + allfit_con[i];
    }
    double* ftemp = new double[n]; // fit of current tree
    // initialization ended

    NumericVector sigma_post(nd);
    NumericVector msd_post(nd);
    NumericVector bsd_post(nd);
    NumericVector b0_post(nd);
    NumericVector b1_post(nd);
    NumericMatrix m_post(nd,n);
    NumericMatrix yhat_post(nd,n);
    NumericMatrix b_post(nd,n);
    IntegerMatrix varcnt_con(nd,p_con);
    IntegerMatrix varcnt_mod(nd,p_mod);
    IntegerMatrix glabel_post(nd,n);
    IntegerMatrix varcnt_cure(nd, p_cure);
    NumericMatrix cure_post(nd,n);
    NumericVector cure_intc_post(nd);
    NumericMatrix pred_cureprob_1(nd,n);
    NumericMatrix pred_cureprob_0(nd,n);

    int save_tree_precision = 32;
    // save stuff to tree file;
    treef_con << std::setprecision(save_tree_precision) << xi_con <<endl;
    treef_con << ntree_con << endl;
    treef_con << di_con.p << endl;
    treef_con << nd << endl;

    treef_mod << std::setprecision(save_tree_precision) << xi_mod << endl;
    treef_mod << ntree_mod << endl;
    treef_mod << di_mod.p << endl;
    treef_mod << nd << endl;
                       
    treef_cure << std::setprecision(save_tree_precision) << xi_cure <<endl;
    treef_cure << ntree_cure << endl;
    treef_cure << di_cure.p << endl;
    treef_cure << nd << endl;

    /* ------------ MCMC ----------------*/

    Rcout << "\n====================================\nBeginning MCMC:\n====================================\n";
    time_t tp;
    int time1 = time(&tp);

    size_t save_ctr = 0;
    bool verbose_itr = false;

    double* weight = new double[n];
    double* weight_heter = new double[n];

    logger.setLevel(0);

    bool printTrees = false;
    
    double u = 0.0; int nuncured = 0; int dicount = 0;

    for (size_t iIter=0; iIter<(nd*thin+burn); iIter++) {
        //verbose_itr = true;
        verbose_itr = false;//iIter>=burn;
        printTrees = false;//iIter>=burn;

        if (verbose_sigma) {
            if (iIter%printevery==0) {
                Rcout << "iteration: " << iIter << " sigma: " << sigma << endl;
            }
        }

        logger.setLevel(verbose_itr);

        logger.log("=========================================");
        sprintf(logBuff, "MCMC iteration: %d of %d Start", iIter+1, nd*thin+burn);
        logger.log(logBuff);
        sprintf(logBuff, "sigma %f", sigma);
        logger.log(logBuff);
        logger.log("==========================================");

        if (verbose_itr) {
            logger.getVectorHead(y, logBuff);
            Rcout << "            y: " << logBuff << "\n";

            logger.getVectorHead(allfit, logBuff);
            Rcout << "Current Fit : " << logBuff << "\n";

            logger.getVectorHead(allfit_con, logBuff);
            Rcout << "allfit_con : " << logBuff << "\n";

            logger.getVectorHead(allfit_mod, logBuff);
            Rcout << "allfit_mod : " << logBuff << "\n";
            
            logger.getVectorHead(glabel, logBuff);
            Rcout << "cured label : " << logBuff << "\n";
            
            logger.getVectorHead(g, logBuff);
            Rcout << "cured latent : " << logBuff << "\n";
            
            logger.getVectorHead(allfit_cure, logBuff);
            Rcout << "allfit_cure : " << logBuff << "\n";
        }
            
            
        logger.log("===============================================");
        logger.log("-- updating latent variables for cured model --");
        logger.log("===============================================");
        
        for (size_t k=0; k<n; k++) {
            if (glabel[k]>0) {
                g[k] = rtnormlo1(allfit_cure[k],0);
            } else {
                g[k] = rtnormhi1(allfit_cure[k],0);
            }
        }
    
        for (size_t k=0; k<n; ++k) {
            weight[k] = w[k]*mscale*mscale/(sigma*sigma);
        }
            
        logger.log("=====================================");
        logger.log("-- Tree Processing for cured model --");
        logger.log("=====================================");
            
        // draw trees for cured trees
        for (size_t iTreecure=0; iTreecure<ntree_cure; iTreecure++) {
            
            logger.log("=======================================");
            sprintf(logBuff, "Updating Cure Tree: %d of %d", iTreecure+1, ntree_cure);
            logger.log(logBuff);
            logger.log("=======================================");
            logger.startContext();
            
            logger.log("Attempting to Print Tree pre Update \n");
            if (verbose_itr && printTrees) {
                t_cure[iTreecure].pr(xi_cure);
                Rcout << "\n\n";
            }
            
            fit(t_cure[iTreecure], xi_cure, di_cure, ftemp);// fit is actually getting mu of each observation and stored in ftemp
            
            logger.log("Attempting to Print Tree Post first call to fit \n");
            if (verbose_itr && printTrees) {
                t_cure[iTreecure].pr(xi_cure);
                Rcout << "\n\n";
            }
            
            for (size_t k=0; k<n; k++) {
                if (ftemp[k] != ftemp[k]) {
                    Rcout << "cured tree " << iTreecure << " obs " << k << " " << endl;
                    Rcout << t_cure[iTreecure] << endl;
                    stop("nan in ftemp");
                }

                allfit_cure[k] = allfit_cure[k] - ftemp[k];

                r_cure[k] = g[k]-allfit_cure[k];

                if (r_cure[k] != r_cure[k]) {
                    Rcout << (g[k]-allfit_cure[k]) << endl;
                    Rcout << r_cure[k] << endl;
                    stop("NaN in resid");
                }
            }
            
            if (verbose_itr && printTrees) {
                logger.getVectorHead(w, logBuff);
                Rcout << "\n weight: " << logBuff << "\n\n";
            }
            
            u = gen.uniform();
            if (u<pi_cure.pbd) {
                logger.log("Starting Birth/Death Processing");
                logger.startContext();
                bd(t_cure[iTreecure], xi_cure, di_cure, w, pi_cure, gen, logger, nv_cure);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post bd \n");
                if (verbose_itr && printTrees) {
                    t_cure[iTreecure].pr(xi_cure);
                    Rcout << "\n";
                }
            } else if(u < (pi_cure.pswap + pi_cure.pbd)) {
                logger.log("Starting SwapRule Processing");
                logger.startContext();
                swaprule(t_cure[iTreecure], xi_cure, di_cure, w, pi_cure, gen, logger);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post SwapRule \n");
                if (verbose_itr && printTrees) {
                    t_cure[iTreecure].pr(xi_cure);
                    Rcout << "\n";
                }
            } else {
                logger.log("Starting ChangeRule Processing");
                logger.startContext();
                changerule(t_cure[iTreecure], xi_cure, di_cure, w, pi_cure, gen, logger, nv_cure);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post ChangeRule \n");
                if (verbose_itr && printTrees) {
                    t_cure[iTreecure].pr(xi_cure);
                    Rcout << "\n";
                }
            }
            
            if (verbose_itr && printTrees) {
                logger.log("Printing Current Status of Fit");

                logger.getVectorHead(a_, logBuff);
                Rcout << "\n           a : " << logBuff << "\n";
                
                logger.getVectorHead(glabel, logBuff);
                Rcout << "\n cured label : " << logBuff << "\n";
                
                logger.getVectorHead(g, logBuff);
                Rcout << "\n cured latent : " << logBuff << "\n";

                logger.getVectorHead(allfit_cure, logBuff);
                Rcout << "Current Fit : " << logBuff << "\n";

                logger.getVectorHead(r_cure, logBuff);
                Rcout << "      r_cure : " << logBuff << "\n\n";
            }

            logger.log("Strarting to draw mu");
            logger.startContext();
            
            drmu(t_cure[iTreecure], xi_cure, di_cure, pi_cure, w, gen);

            logger.stopContext();

            logger.log("Attempting to Print Tree Post drmu \n") ;
            if (verbose_itr && printTrees) {
                t_cure[iTreecure].pr(xi_cure);
                Rcout << "\n";
            }

            fit(t_cure[iTreecure], xi_cure, di_cure, ftemp);

            for (size_t k=0; k<n; k++) {
                allfit_cure[k] += ftemp[k];
            }

            logger.log("Attempting to Print Tree Post second call to fit \n");

            if (verbose_itr && printTrees) {
                t_cure[iTreecure].pr(xi_cure);
                Rcout << "\n";
            }
            logger.stopContext();
        } // end tree loop for the cured model
            
        // pick the uncured subjects and re-arrange used for AFT model through the weights
        nuncured = 0;
        for (size_t k=0; k<ntrt; k++) {
            if (glabel[k]>0) {
                nuncured += 1;
                weight[k] = w[k]*mscale*mscale/(sigma*sigma);
                weight_heter[k] = w[k]*bscale1*bscale1/(sigma*sigma);
            } else {
                weight[k] = w[k]*mscale*mscale/(sigma*sigma);//0.;
                weight_heter[k] = w[k]*bscale1*bscale1/(sigma*sigma);//0.;
            }
        }
            
        for (size_t k = ntrt; k<n; k++) {
            if (glabel[k]>0) {
                nuncured += 1;
                weight[k] = w[k]*mscale*mscale/(sigma*sigma);
                weight_heter[k] = w[k]*bscale0*bscale0/(sigma*sigma);
            } else {
                weight[k] = w[k]*mscale*mscale/(sigma*sigma);//0.;
                weight_heter[k] = w[k]*bscale0*bscale0/(sigma*sigma);//0.;
            }
        }
        
        
        double* r_con_temp = new double[nuncured];
        double* x_con_temp = new double[nuncured*p_con];
        double* r_mod_temp = new double[nuncured];
        double* x_mod_temp = new double[nuncured*p_mod];
        double* weight_temp = new double[nuncured];
        double* weight_heter_temp = new double[nuncured];
        dicount = 0;
        for (size_t k=0; k<n; ++k) {
            if (glabel[k] > 0) {
                weight_temp[dicount] = weight[k];
                weight_heter_temp[dicount] = weight_heter[k];
                for (size_t p=0; p<p_con; ++p) {
                    x_con_temp[dicount*p_con + p] = x_con[k*p_con + p];
                }
                for (size_t p=0; p<p_mod; ++p) {
                    x_mod_temp[dicount*p_mod + p] = x_mod[k*p_mod + p];
                }
                dicount += 1;
            }
        }
        
    

        logger.log("==========================");
        logger.log("-- Tree Processing --");
        logger.log("==========================");

        // draw trees for m(x);
        for (size_t iTreecon=0; iTreecon<ntree_con; iTreecon++) {

            logger.log("===========================");
            sprintf(logBuff, "Updating Control Tree: %d of %d", iTreecon+1, ntree_con);
            logger.log(logBuff);
            logger.log("===========================");
            logger.startContext();

            logger.log("Attempting to Print Tree pre Update \n");
            if (verbose_itr && printTrees) {
                t_con[iTreecon].pr(xi_con);
                Rcout << "\n\n";
            }

            
            fit(t_con[iTreecon], xi_con, di_con, ftemp);

            logger.log("Attempting to Print Tree Post first call to fit \n");
            if (verbose_itr && printTrees) {
                t_con[iTreecon].pr(xi_con);
                Rcout << "\n\n";
            }

            for (size_t k=0; k<n; k++) {
                if (ftemp[k] != ftemp[k]) {
                    Rcout << "control tree " << iTreecon << " obs " << k << " " << endl;
                    Rcout << t_con[iTreecon] << endl;
                    stop("nan in ftemp");
                }
            }
            
            
            dicount = 0;
            for (size_t k=0; k<n; k++) {
                allfit[k] = allfit[k] - mscale*ftemp[k];
                allfit_con[k] = allfit_con[k] - mscale*ftemp[k];
                r_con[k] = (y[k]-allfit[k])/mscale;
                if (r_con[k] != r_con[k]) {
                    Rcout << (y[k]-allfit[k]) << endl;
                    Rcout << mscale << endl;
                    Rcout << r_con[k] << endl;
                    stop("NaN in resid");
                }
                if (glabel[k]>0) {
                    r_con_temp[dicount] = r_con[k];
                    dicount += 1;
                }
            }


            if (verbose_itr && printTrees) {
                logger.getVectorHead(weight_temp, logBuff);
                Rcout << "\n weight: " << logBuff << "\n\n";
            }
            
            di_con.N = nuncured;
            di_con.x = x_con_temp;
            di_con.y = r_con_temp;
            
            
            u = gen.uniform();
            if (u<pi_con.pbd) {
                logger.log("Starting Birth/Death Processing");
                logger.startContext();
                bd(t_con[iTreecon], xi_con, di_con, weight_temp, pi_con, gen, logger, nv_con);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post bd \n");
                if (verbose_itr && printTrees) {
                    t_con[iTreecon].pr(xi_con);
                    Rcout << "\n";
                }
            } else if(u < (pi_con.pswap + pi_con.pbd)) {
                logger.log("Starting SwapRule Processing");
                logger.startContext();
                swaprule(t_con[iTreecon], xi_con, di_con, weight_temp, pi_con, gen, logger);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post SwapRule \n");
                if (verbose_itr && printTrees) {
                    t_con[iTreecon].pr(xi_con);
                    Rcout << "\n";
                }
            } else {
                logger.log("Starting ChangeRule Processing");
                logger.startContext();
                changerule(t_con[iTreecon], xi_con, di_con, weight_temp, pi_con, gen, logger, nv_con);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post ChangeRule \n");
                if (verbose_itr && printTrees) {
                    t_con[iTreecon].pr(xi_con);
                    Rcout << "\n";
                }
            }
            

            if (verbose_itr && printTrees) {
                logger.log("Printing Current Status of Fit");

                logger.getVectorHead(a_, logBuff);
                Rcout << "\n           a : " << logBuff << "\n";
                
                logger.getVectorHead(y, logBuff);
                Rcout << "\n           y : " << logBuff << "\n";

                logger.getVectorHead(allfit, logBuff);
                Rcout << "Current Fit : " << logBuff << "\n";

                logger.getVectorHead(r_con, logBuff);
                Rcout << "      r_con : " << logBuff << "\n\n";
            }

            logger.log("Strarting to draw mu");
            logger.startContext();

            drmu(t_con[iTreecon], xi_con, di_con, pi_con, weight_temp, gen);

            logger.stopContext();

            logger.log("Attempting to Print Tree Post drmu \n") ;
            if (verbose_itr && printTrees) {
                t_con[iTreecon].pr(xi_con);
                Rcout << "\n";
            }

            di_con.N = n;
            di_con.x = &x_con[0];
            fit(t_con[iTreecon], xi_con, di_con, ftemp);
            
            
            
           
            for (size_t k=0; k<n; k++) {
                allfit[k] += mscale*ftemp[k];
                allfit_con[k] += mscale*ftemp[k];
            }

            logger.log("Attempting to Print Tree Post second call to fit \n");

            if (verbose_itr && printTrees) {
                t_con[iTreecon].pr(xi_con);
                Rcout << "\n";
            }
            logger.stopContext();
        } // end tree loop


        // Next, for b(x) trees

        for (size_t iTreeMod=0; iTreeMod<ntree_mod; iTreeMod++) {
            logger.log("============================");
            sprintf(logBuff,"Updating Moderate Tree: %d of %d", iTreeMod+1, ntree_mod);
            logger.log(logBuff);
            logger.log("============================");
            logger.startContext();

            logger.log("Attempting to Print Tree Pre Update \n");
            if(verbose_itr && printTrees){
                    t_mod[iTreeMod].pr(xi_mod);
                    Rcout << "\n";
                  }
            

            fit(t_mod[iTreeMod],
                      xi_mod,
                      di_mod,
                      ftemp);

            logger.log("Attempting to Print Tree Post first call to fit");
            if(verbose_itr && printTrees){
                    t_mod[iTreeMod].pr(xi_mod);
                    Rcout << "\n";
                  }

            for (size_t k=0; k<n; k++) {
                if (ftemp[k] != ftemp[k]) {
                    Rcout << "moderate tree " << iTreeMod << " obs " << k << " " << endl;
                    Rcout << t_mod[iTreeMod] << endl;
                    stop("nan in ftemp");
                }
            }
            
            dicount = 0;
            for (size_t k=0; k<n; k++) {
                double bscale = (k<ntrt) ? bscale1 : bscale0;
                allfit[k] = allfit[k] - bscale*ftemp[k];
                allfit_mod[k] = allfit_mod[k] - bscale*ftemp[k];
                r_mod[k] = (y[k] - allfit[k])/bscale;
                if (glabel[k]>0) {
                    r_mod_temp[dicount] = r_mod[k];
                    dicount += 1;
                }
            }
            
            di_mod.N = nuncured;
            di_mod.x = x_mod_temp;
            di_mod.y = r_mod_temp;
            
            u = gen.uniform();
            if (u<pi_mod.pbd) {
                logger.log("Starting Birth/Death Processing");
                logger.startContext();
                bd(t_mod[iTreeMod], xi_mod, di_mod, weight_heter_temp, pi_mod, gen, logger, nv_mod);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post bd \n");
                if (verbose_itr && printTrees) {
                    t_mod[iTreeMod].pr(xi_mod);
                    Rcout << "\n";
                }
            } else if (u < (pi_mod.pbd + pi_mod.pswap)) {
                logger.log("Starting SwapRule Processing");
                logger.startContext();
                swaprule(t_mod[iTreeMod], xi_mod, di_mod, weight_heter_temp, pi_mod, gen, logger);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post SwapRule \n");
                if (verbose_itr && printTrees) {
                    t_mod[iTreeMod].pr(xi_mod);
                    Rcout << "\n";
                }
            } else {
                logger.log("Starting ChangeRule Processing");
                logger.startContext();
                changerule(t_mod[iTreeMod], xi_mod, di_mod, weight_heter_temp, pi_mod, gen, logger, nv_mod);
                logger.stopContext();

                logger.log("Attempting to Print Tree Post ChangeRule \n");
                if (verbose_itr && printTrees) {
                    t_mod[iTreeMod].pr(xi_mod);
                    Rcout << "\n";
                }
            }
            
            
            if (verbose_itr && printTrees) {
                logger.log("Printing Status of Fit");

                logger.getVectorHead(a_, logBuff);
                Rcout << "\n       a : " << logBuff << "\n";
                
                logger.getVectorHead(y, logBuff);
                Rcout << "         y : " << logBuff << "\n";

                logger.getVectorHead(allfit, logBuff);
                Rcout << " Fit - Tree : " << logBuff << "\n";
                logger.getVectorHead(r_mod, logBuff);
                Rcout << "    r_mod : " << logBuff << "\n\n";
                Rcout << "mscale: " << mscale << "\n";

                Rcout << "bscale0: " << bscale0 << "\n";

                Rcout << "bscale1: " << bscale1 << "\n\n";
            }
            logger.log("Starting to draw mu");
            logger.startContext();
            drmu(t_mod[iTreeMod], xi_mod, di_mod, pi_mod, weight_heter_temp, gen);
            logger.stopContext();

            logger.log("Attempting to Print Tree Post drmuheter \n");
            if (verbose_itr && printTrees) {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }

            di_mod.N = n;
            di_mod.x = &x_mod[0];
            
            fit(t_mod[iTreeMod], xi_mod, di_mod, ftemp);

            for (size_t k=0; k<ntrt; k++) {
                allfit[k] += bscale1*ftemp[k];
                allfit_mod[k] += bscale1*ftemp[k];
            }
            
            for (size_t k=ntrt; k<n; k++) {
                allfit[k] += bscale0*ftemp[k];
                allfit_mod[k] += bscale0*ftemp[k];
            }

            logger.log("Attempting to Print Tree Post second call to fit");

            if (verbose_itr && printTrees) {
                t_mod[iTreeMod].pr(xi_mod);
                Rcout << "\n";
            }
            logger.stopContext();
        } // end tree loop

        logger.setLevel(verbose_itr);

        logger.log("============================");
        logger.log("-- MCMC iteration cleanup --");
        logger.log("============================");

        
        // -----------------
        logger.log("Draw Sigma");
        // ------------------
        double rss = 0.0;
        double restemp = 0.0;
        for (size_t k=0; k<n; k++) {
            if (glabel[k]>0) {
                restemp = y[k] - allfit[k];
                rss += w[k]*restemp*restemp;
            }
        }
        sigma = sqrt((nu*kappa + rss)/gen.chisq(nu+nuncured));
        pi_con.sigma = sigma/fabs(mscale);
        pi_mod.sigma = sigma;
            
        // for cencored subjects
        if (num_censored > 0) {
            for (size_t k=0; k<n; k++) {
                if ((delta[k]==0) & (glabel[k]>0)) {
                    y[k] = rtnormlo(allfit[k], sigma, yobs[k]);
                    //sprintf(logBuff, "%d of %d censored subjects imputed, yobs : %f, yimp : %f", k, num_censored, yobs[k], y[k]);
                    //Rcout << logBuff << endl;
                }
            }
        }
            
        // update the cured group labels
        for (size_t k=0; k<n; k++) {
            if (delta[k]==0) {
                double p1c = R::pnorm(allfit_cure[k], 0., 1., true, false);
                double p2s = R::pnorm(yobs[k] - allfit[k], 0., sigma, false, false);
                double p_post = (p1c*p2s)/(1-p1c+p1c*p2s);
                glabel[k] = gen.binom(1,p_post);//(int) R::rbinom(1,p_post);
            }
        }
            
            
            
        // store posterior draws
        if (((iIter>=burn) & (iIter % thin==0))) {
            for (size_t j=0; j<ntree_con; j++) treef_con << std::setprecision(save_tree_precision) << t_con[j] << endl; // save trees
            for (size_t j=0; j<ntree_mod; j++) treef_mod << std::setprecision(save_tree_precision) << t_mod[j] << endl;
            for (size_t j=0; j<ntree_cure; j++) treef_cure << std::setprecision(save_tree_precision) << t_cure[j] << endl;

            msd_post(save_ctr) = mscale;
            bsd_post(save_ctr) = bscale1 - bscale0;
            b0_post(save_ctr) = bscale0;
            b1_post(save_ctr) = bscale1;

            sigma_post(save_ctr) = sigma;
            for (size_t k=0; k<n; k++) {
                m_post(save_ctr,k) = allfit_con[k];
                yhat_post(save_ctr,k) = allfit[k];
            }
            for (size_t k=0; k<n; k++) {
                double bscale = (k<ntrt) ? bscale1 : bscale0;
                b_post(save_ctr, k) = (bscale1 - bscale0)*allfit_mod[k]/bscale;
            }
            
            for (size_t k=0; k<p_con; k++) {
                varcnt_con(save_ctr,k) = nv_con[k];
            }
            
            for (size_t k=0; k<p_mod; k++) {
                varcnt_mod(save_ctr,k) = nv_mod[k];
            }
            
            for (size_t k=0; k<p_cure; k++) {
                varcnt_cure(save_ctr,k) = nv_cure[k];
            }
            
            for (size_t k=0; k<n; k++) {
                glabel_post(save_ctr,k) = glabel[k];
                cure_post(save_ctr,k) = allfit_cure[k];
            }
            
            // do the prediction
            predict(t_cure, xi_cure, di_pred_cure_0, fitpred_cure_0);
            predict(t_cure, xi_cure, di_pred_cure_1, fitpred_cure_1);
            for (size_t k=0; k<n; k++) {
                pred_cureprob_0(save_ctr,k) = R::pnorm(fitpred_cure_0[k], 0., 1., true, false);
                pred_cureprob_1(save_ctr,k) = R::pnorm(fitpred_cure_1[k], 0., 1., true, false);
            }

            
            save_ctr += 1;
            //treef_con << std::setprecision(save_tree_precision) << save_ctr << endl;
            //treef_mod << std::setprecision(save_tree_precision) << save_ctr << endl;
            
        }


        logger.log("===================================");
        sprintf(logBuff, "MCMC iteration: %d of %d End", iIter+1, nd*thin+burn);
        logger.log(logBuff);
        sprintf(logBuff, "sigma %f, mscale %f, bscale0 %f, bscale1 %f", sigma, mscale, bscale0, bscale1);
        logger.log(logBuff);
        logger.log("===================================");

        if (verbose_itr) {
            logger.getVectorHead(y, logBuff);
            Rcout << "      y : " << logBuff << "\n";
            
            logger.getVectorHead(allfit, logBuff);
            Rcout << " Current Fit: " << logBuff << "\n";

            logger.getVectorHead(allfit_con, logBuff);
            Rcout << "allfit_con: " << logBuff << "\n";

            logger.getVectorHead(allfit_mod, logBuff);
            Rcout << "allfit_mod: " << logBuff << "\n";
        }
        
        delete[] r_con_temp;
        delete[] r_mod_temp;
        delete[] x_con_temp;
        delete[] x_mod_temp;
        delete[] weight_temp;
        delete[] weight_heter_temp;


    } // end MCMC loop;

    int time2 = time(&tp);
    Rcout << "\n=========================\n MCMC Complete \n=========================\n";
    Rcout << "time for loop: " << time2 - time1 << endl;

    t_mod.clear(); t_con.clear();t_cure.clear();
    delete [] allfit;
    delete[] allfit_mod;
    delete[] allfit_con;
    delete[] r_mod;
    delete[] r_con;
    delete[] ftemp;
    delete [] fitpred_cure_1;
    delete [] fitpred_cure_0;
    delete [] allfit_cure;
    delete [] r_cure;

    treef_con.close();
    treef_mod.close();
    treef_cure.close();

    return (List::create(_["yhat_post"] = yhat_post, _["m_post"] = m_post, _["b_post"] = b_post, _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post, _["b0"] = b0_post, _["b1"] = b1_post, _["varcnt_con_post"] = varcnt_con, _["varcnt_mod_post"] = varcnt_mod, _["varcnt_cure_post"] = varcnt_cure, _["glabel_post"] = glabel_post, _["curefit_post"] = cure_post, _["pred_cureprob_1_post"] = pred_cureprob_1, _["pred_cureprob_0_post"] = pred_cureprob_0));


}



