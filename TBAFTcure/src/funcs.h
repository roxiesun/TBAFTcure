#ifndef GUARD_funcs_h
#define GUARD_funcs_h

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rand_gen.h"

inline double logsumexp(const double &a, const double &b){
    return a<b? b + log(1.0 + exp(a-b)): a + log(1.0 + exp(b-a));
}

using std::cout;
using std::endl;

// define log2pi
#define LTPI 1.837877066

typedef std::vector<std::vector<int> > lookup_t;

lookup_t make_lookup(Rcpp::IntegerMatrix lookup_table, Rcpp::IntegerVector cx);

void impute_x(int v,
              std::vector<int>& mask,
              int n, xinfo& xi, std::vector<double>& x, std::vector<vector<int> >&x_cat,
              std::vector<int>& cx, std::vector<int>& offsets, std::vector<int>& x_type,
              std::vector<tree>& t, std::vector<double>& y, double& sigma, RNG& rng);

//normal density
double pn(double x, double m, double v);

// discrete distribution draw
int rdisc(double *p, RNG& gen);

// evaluate tree on grid xinfo, and write
void grm(tree&tr, xinfo& xi, std::ostream& os);

// check whether a node has variables it can split on
bool cansplit(tree::tree_p n, xinfo& xi);

// compute prob birth
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots);

// find variables can split on and store in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);

void getinternalvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars);

int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);

//-------------------
// get prob a node grows, 0 if no good vars, a/(1+d)^b else
double pgrow(tree::tree_p n, xinfo&xi, pinfo& pi);

// Calibart
void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U);

void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U);

//-------------------
// suff statistics for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, double* weight, tree::npv& bnv, std::vector<sinfo>& sv);
// counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di);
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv);

// update counts to reflect observation i?
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, int sign);

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, int sign);

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp, size_t>& bnmap, int sign);

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp, size_t>& bnmap, int sign, tree::tree_cp &tbn);
//-------------------

// check min leaf size
bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di);

//-------------------
// get sufficient stat for children (v,c) of node nx in tree x
void getsuffBirth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);

// get sufficient stat for pair of bottom children nl, nr, in tree x
void getsuffDeath(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr);

//---------------------
// log integrated likelihood
double loglike(double n, double sy, double sigma, double tau);

//log pri of Tree
double logPriT(tree::tree_p x, xinfo& xi, pinfo& pi);
// fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv);

void fit(tree& t, xinfo& xi, dinfo& di, double* fv);

void predict(std::vector<tree>& tv, xinfo& xi, dinfo& di, double* fp);

// template function?
// one tree nodes mean?
template<class T>
double fit_i(T i, tree& t, xinfo& xi, dinfo& di)
{
    double *xx;
    double fv = 0.0;
    tree::tree_cp bn;
    xx = di.x + i*di.p; //what for?
    bn = t.bot(xx,xi); // find leaf
    fv = bn->getmu();
    
    return fv;
}

// forest nodes mean?
template<class T>
double fit_i(T i, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
    double *xx;
    double fv = 0.0;
    tree::tree_cp bn;
    xx = di.x + i*di.p;
    for (size_t j=0; j<t.size(); ++j) {
        bn = t[j].bot(xx,xi);
        fv += bn->getmu();
    }
    return fv;
}

template<class T>
double fit_i_mult(T i, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
    double *xx;
    double fv = 1.0;
    tree::tree_cp bn;
    xx = di.x + i*di.p;
    for (size_t j=0; j<t.size(); ++j) {
        bn = t[j].bot(xx,xi);
        fv *= bn->getmu();
    }
    return fv;
}

//---------------------
//partition

void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv);

// draw all bottom node mu's
#ifdef MPIBART

void MPImasterdrmu(tree& t, xinfo& xi, pinfo&pi, RNG& gen, size_t numslaves);
#else
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double* weight, RNG& gen);
void drphi(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
#endif

//---------------------
// write xinfo to screen
void prxi(xinfo& xi);

//---------------
//make xinfo
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
// get min/max for p predictors needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx);
// make xinfo given min/max
void makexinfominmax(size_t p, xinfo&xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx);

// fucntions for updating mixture quatities
void updateLabels(int* labels, double* mixprop, double* locations, double* resid, double sigma, size_t nobs, size_t nmix, RNG& gen);
void updateMixprp(double* mixprop, double* mass, int* labels, int* mixcnts, size_t nmix, double psi1, double psi2, RNG& gen);
void updateLocations(double* locations, double* mixprop, double* resid, int* labels, int* mixcnts, double sigma, double prior_sigsq, size_t nobs, size_t nmix, RNG& gen);
void updateIndivLocations (double* indiv_locations, double* locations, int* labels, size_t nobs, size_t nmix);
bool CheckRule(tree::tree_p n, xinfo& xi, size_t var);

#ifdef MPIBART
//MPI calls
void MPImasterallsuff(tree& x, tree::npv& bnv, std::vector<sinfo>& sv, size_t numslaves);
void MPIslaveallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv);
void MPIslavedrmu(tree& t, xinfo& xi, dinfo& di);
void MPImastergetsuff(tree::tree_cp nl, tree::tree_cp nr, sinfo &sl, sinfo &sr, size_t numslaves);
void MPImastergetsuffvc(tree::tree_cp nx, size_t v, size_t c, xinfo& xi, sinfo& sl, sinfo& sr, size_t numslaves);
void MPImastersendbirth(tree::tree_p nx, size_t v, size_t c, double mul, double mur, size_t numslaves);
void MPImastersenddeath(tree::tree_p nx, double mu, size_t numslaves);
void MPImastersendnobirthdeath(size_t numslaves);
void MPIslaveupdatebirthdeath(tree& x);
void MPIslavegetsuff(tree& x, xinfo& xi, dinfo& di);
void makepred(dinfo dip, xinfo &xi, std::vector<tree> &t, double *ppredmean);
void makeypred(dinfo dip, xinfo &xi, std::vector<tree> &t, double sigma, double *ppredmean);
void makepostpred(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp);
void makepostpred2(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp, double *postp2);
double logcalp(std::vector<double> &theta, dinfo dip, xinfo &xi, std::vector<tree> &t, double sigmae, double sigma, size_t pth, size_t myrank);
#endif

#endif /* funcs_h */
