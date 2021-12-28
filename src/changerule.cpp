#include <iostream>
#include <RcppArmadillo.h>

#include "info.h"
#include "tree.h"
#include "changerule.h"
#include "funcs.h"

using std::cout;
using std::endl;

bool changerule(tree& x, xinfo& xi, dinfo& di, double*phi, pinfo& pi, RNG& gen, Logger logger, std::vector<size_t>& ivcnt)
{
    tree::npv internv; // not bot nodes
    x.getnbots(internv);
    size_t Ninternv = internv.size();
    logger.log("Attempting ChangeRule");
    if (Ninternv == 0) {
        logger.log("Rejecting ChangeRule");
        return false;
    }
    
    size_t ni = floor(gen.uniform()*Ninternv);//randomly choose a not bot nodes to change its rule;
    tree::tree_p nx = internv[ni];
    
    // find a good var for the chosen INTERNAL node nx
    std::vector<size_t> goodintvars;
    getinternalvars(nx, xi, goodintvars);
    size_t vi = floor(gen.uniform()*goodintvars.size());
    size_t v = goodintvars[vi];
    
    // draw cutpoint uniformly for node nx with var v
    int L,U;
    L=0; U=xi[v].size()-1;
    //nx->region(v, &L, &U); not true, because our internal node may have left/right children that also uses v,instead, we use
    getpertLU(nx, v, xi, &L, &U);
    size_t c = L + floor(gen.uniform()*(U-L+1));
    
    // prepare for MH
    // note that the proposal distribution is identical for the changerule move??, so we only need to consider the logIlike,
    std::vector<sinfo> oldsv;
    tree::npv bnv;
    allsuff(x, xi, di, phi, bnv, oldsv);
    double logILold = 0.0;
    for (tree::npv::size_type i=0 ; i!=bnv.size(); i++) {
        logILold += loglike(oldsv[i].n, oldsv[i].sy, pi.sigma, pi.tau);
    }
    double logPriold = logPriT(nx, xi, pi);
    
    // save current var and cutoff
    size_t oldv = nx->getv();
    size_t oldc = nx->getc();
    
    // change to new rule
    tree::tree_p np = x.getptr(nx->nid());
    np->setv(v);
    np->setc(c);
    // because I doubt nx is just a copy of orginal tree node and thus setting v&c for nx might not work.
    //nx->setv(v);
    //nx->setc(c);
    std::vector<sinfo> newsv;
    tree::npv bnvnew;
    allsuff(x, xi, di, phi, bnvnew, newsv);
    double logInew = 0.0;
    for (tree::npv::size_type i=0; i!=bnvnew.size(); i++) {
        logInew += loglike(newsv[i].n, newsv[i].sy, pi.sigma, pi.tau);
    }
    double logPrinew = logPriT(np, xi, pi);
    
    double alpha = std::min(1.0, exp(logPrinew + logInew - logPriold - logILold));
    
    if (gen.uniform() < alpha) {
        logger.log("Accepting ChangeRule");
        ivcnt[oldv] -= 1;
        ivcnt[v] += 1;
        return true;
    } else {
        np->setv(oldv);
        np->setc(oldc);
        logger.log("Rejecting ChangeRule");
        return  false;
    }
}
