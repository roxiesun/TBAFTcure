#include <iostream>
#include <RcppArmadillo.h>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funcs.h"

using std::cout;
using std::endl;

/*Note: going from state x to state y
 if return true: birth/death accept*/

bool bd(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen, Logger logger, std::vector<size_t>& ivcnt)
{
    tree::npv goodbots; // bot nodes that can split
    double PBx = getpb(x, xi, pi, goodbots);
    
    if (gen.uniform() < PBx) {
        logger.log("Attempting Birth");
        
        // draw proposal
        
        // uniformly draw bottom node: choose node index from goodbots
        size_t ni = floor(gen.uniform()*goodbots.size()); // rounddown
        tree::tree_p nx = goodbots[ni];
        
        // draw variable v, uniformly
        std::vector<size_t> goodvars;
        getgoodvars(nx, xi, goodvars); // get variable this node can split on
        size_t vi = floor(gen.uniform()*goodvars.size());
        size_t v = goodvars[vi];
        
        // draw cutpoint, uniformly
        int L,U;
        L=0; U=xi[v].size()-1;
        nx->region(v, &L, &U);
        size_t c = L + floor(gen.uniform()*(U-L+1));
        // U-L+1 is the number of available split points
        
        //-------------------
        // prepare for Metropolis hastings
        double Pbotx = 1.0/goodbots.size(); // proposal dist/probability of choosing nx;
        size_t dnx = nx->depth();
        double PGnx = pi.alpha/pow(1.0+dnx, pi.beta); // prior prob of growing at nx;
        
        double PGly,PGry; // prior probs of growing at new children
        if (goodvars.size()>1) {
            PGly = pi.alpha/pow(1.0+dnx+1.0, pi.beta);
            PGry = PGly;
        } else { // have only one v to work with
            if ((int)(c-1)<L) { // v exhausted in new left child l, new upper limit would be c-1
                PGly = 0.0;
            } else {
                PGly = pi.alpha/pow(1.0+dnx+1.0, pi.beta);
            }
            if (U<(int)(c+1)) { // v exhausted in new right child r, new lower limit would be c+1
                PGry = 0.0;
            } else {
                PGry = pi.alpha/pow(1.0+dnx+1.0, pi.beta);
            }
        }
        
        double PDy; // prob of proposing death at y;
        if (goodbots.size()>1) { // can birth at y because splittable nodes left
            PDy = 1.0 - pi.pb;
        } else { //nx is the only splittable node
            if ((PGry==0) && (PGly==0)) { //cannot birth at y
                PDy = 1.0;
            } else { // y can birth can either l or r
                PDy = 1.0 - pi.pb;
            }
        }
        
        double Pnogy; // death prob of choosing the nog node at y
        size_t nnogs = x.nnogs();
        tree::tree_cp nxp = nx->getp();
        if (nxp==0) {
            Pnogy = 1.0;
        } else {
            if (nxp->isnog()) { // is parent is a nog, nnogs same at state x and y
                Pnogy = 1.0/nnogs;
            } else {
                Pnogy = 1.0/(nnogs+1.0);
            }
        }
        
        //-----------
        //compute suff stats
        sinfo sl,sr;
        getsuffBirth(x, nx, v, c, xi, di, phi, sl, sr);
        
        //---------
        //compute alpha
        double alpha=0.0, alpha1=0.0, alpha2=0.0;
        double lill=0.0, lilr=0.0, lilt=0.0;
        if ((sl.n0>4) && (sr.n0>4)) {
            lill = loglike(sl.n, sl.sy, pi.sigma, pi.tau);
            lilr = loglike(sr.n, sr.sy, pi.sigma, pi.tau);
            lilt = loglike(sl.n+sr.n, sl.sy+sr.sy, pi.sigma, pi.tau);
            
            alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx); // prior*proposal prob is this odd?, yes, after takeing exp this produces the difference in likelihood
            alpha2 = alpha1*exp(lill+lilr-lilt);//prior*proposal*likelihood
            alpha =std::min(1.0, alpha2);
        } else {
            alpha = 0.0;
        }
        
        //--------------
        // finally MH
        double a,b,s2,yb;
        double mul,mur;
        
        if (gen.uniform()<alpha) { //accept birth
            logger.log("Accepting Birth");
            // draw mul
            a = 1.0/(pi.tau*pi.tau); //1/tau^2
            s2 = pi.sigma*pi.sigma; //sigma^2
            //left mean
            yb = sl.sy/sl.n;
            b = sl.n/s2;
            mul = b*yb/(a+b) + gen.normal()/sqrt(a+b);
            //draw mur
            yb = sr.sy/sr.n;
            b = sr.n/s2;
            mur = b*yb/(a+b) + gen.normal()/sqrt(a+b);
            //do birth
            x.birth(nx->nid(), v, c, mul, mur);
            ivcnt[v] += 1;
            return true;
        } else {
            logger.log("Rejecting Birth");
            return false;
        }
    } else { // if not do birth, do death
        logger.log("Attempting Death:");
        
        //draw proposal
        // draw a nog, any nog is possible
        tree::npv nognodes;
        x.getnogs(nognodes);
        size_t ni = floor(gen.uniform()*nognodes.size());
        tree::tree_p nx = nognodes[ni];
        
        // prepare for MH ratio
        
        double PGny; //prob the nog node grows;
        size_t dny = nx->depth();
        PGny = pi.alpha/pow(1.0+dny, pi.beta);
        
        double PGlx = pgrow(nx->getl(), xi, pi);
        double PGrx = pgrow(nx->getr(), xi, pi);
        
        double PBy; // prob of birth move at y
        if (!(nx->p)) { // is the nog nx the top node
            PBy = 1.0;
        } else {
            PBy = pi.pb;
        }
        
        double Pboty; // prob of choosing the nog as bot to split on??
        int ngood = goodbots.size();
        if (cansplit(nx->getl(), xi)) --ngood; // if can split at left child, loose this one;
        if (cansplit(nx->getr(), xi)) --ngood;
        ++ngood; // if we can split can nx
        Pboty = 1.0/ngood;
        
        double PDx = 1.0-PBx; // prob pf death step at x(state);
        double Pnogx = 1.0/nognodes.size();
        
        // suff stats
        sinfo sl,sr;
#ifdef MPIBART
        MPImastergetsuff(x, nx->getl(),nx->getr(),sl,sr,numslaves);
#else
        getsuffDeath(x, nx->getl(), nx->getr(), xi, di, phi, sl, sr);
#endif
        
        //--------------
        //compute alpha
        
        double lill = loglike(sl.n, sl.sy, pi.sigma, pi.tau);
        double lilr = loglike(sr.n, sr.sy, pi.sigma, pi.tau);
        double lilt = loglike(sl.n+sr.n, sl.sy+sr.sy, pi.sigma, pi.tau);
        
        double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
        double alpha2 = alpha1*exp(lilt - lill - lilr);
        double alpha = std::min(1.0, alpha2);
        
        // finally MH
        double a,b,s2,yb;
        double mu;
        double n;
        
        if (gen.uniform()<alpha) { //accept death
            logger.log("Accepting Death");
            //draw mu for nog(which will be bot)
            n = sl.n+sr.n;
            a = 1.0/(pi.tau*pi.tau);
            s2 = pi.sigma*pi.sigma;
            yb = (sl.sy+sr.sy)/n;
            b = n/s2;
            mu = b*yb/(a+b) + gen.normal()/sqrt(a+b);
            // do death;
            ivcnt[nx->getv()] -= 1;
            x.death(nx->nid(), mu);
#ifdef MPIBART
            MPImastersenddeath(nx,mu,numslaves);
#endif
            return true;
        } else {
            logger.log("Rejecting Death");
#ifdef MPIBART
            MPImastersendnobirthdeath(numslaves);
#endif
            return false;
        }
        
    }
}
