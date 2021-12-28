#ifndef GUARD_info_h
#define GUARD_info_h

#include<vector>

// data info, define constants
class dinfo {
public:
    dinfo() {p = 0; N = 0; x = 0; y = 0;} // constructor?
    size_t p;   // dim of covariates
    size_t N;   // number of observations
    double *x;  /* so that j_th var of the i_th subject is *(x + p*i + j), saved as a vector? */
    double *y; // y_i is *(y+i) or y[i]
};

// prior and mcmc information
class pinfo {
public:
    pinfo() {pbd = 0.5; pb = 0.25; pchange = 0.4; pswap = 0.1; alpha = 0.95; beta = 0.5; tau = 1.0; sigma = 1.0; minperbot = 5;
        taccept = 0; tavgd = 0; tmaxd = 0;}
    //mcmc
    double pbd;  // prob of birth/death
    double pb;   // prob of birth
    // But what about swap and change?
    double pchange;
    double pswap;
    
    // prior info
    double alpha; // base hyper
    double beta;  // power hyper?
    double tau;   // true homogenous treatment effects? or what hyper? node variance hyperparameter?
    
    // sigma
    double sigma; //error variance
    
    //tree min obs
    size_t minperbot; // what for?
    
    // mcmc settings
    std::vector< std::vector<double> > cvpm;  // change of variable probability matrix
    unsigned int taccept;  // acceptance count of tree proposals
    unsigned int tproposal; // number of tree proposals, to check the acceptance rate?
    unsigned int tavgd; // average tree depth
    unsigned int tmaxd; // maximum tree depth
};

// sufficient statistics for 1 node
class sinfo {
public:
    sinfo() {n0 = 0.0; n = 0; sy = 0.0;} // constructor
    double n0;  // unweighted sample counts.
    double n;
    double sy;  // sum y of this node?
};

#endif /* info_h */
