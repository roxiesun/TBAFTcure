// [[Rcpp::depends(RcppParallel)]]
#include <cmath>
#include "funcs.h"
#include <map>
#include <RcppParallel.h>
#ifdef MPIBART
#include "mpi.h"
#endif

using Rcpp::Rcout;
using namespace RcppParallel;
using Rcpp::stop;

// normal density
double pn(double x, double m, double v)
{
    double dif = x-m;
    return exp(-0.5*dif*dif/v)/sqrt(2*PI*v);
}

// discrete distribution draw
int rdisc(double *p, RNG& gen)
{
    double sum;
    double u = gen.uniform();
    
    int i = 0;
    sum = p[0];
    while (sum < u) {
        i += 1;
        sum += p[i];
    }
    return i;
}

// evaluate tree on grid xinfo, and write
void grm(tree&tr, xinfo& xi, std::ostream& os)
{
    size_t p = xi.size();
    //check the dim of xi
    if (p!=2) {
        Rcout << "error in grm: p!=2\n";
        return;
    }
    size_t n1 = xi[0].size();
    size_t n2 = xi[1].size();
    tree::tree_cp bp;
    double *x = new double[2]; // array of two elements
    for (size_t i=0; i!=n1; i++) {
        for (size_t j=0; j!=n2; j++) {
            x[0] = xi[0][i];
            x[1] = xi[1][j];
            bp = tr.bot(x, xi);
            os << x[0] << " " << x[1] << " " << bp->getmu() << " " << bp->nid() << endl;
        }
    }
    delete[] x;
}

// check whether a node has variables it can split on
bool cansplit(tree::tree_p n, xinfo& xi)
{
    int L, U;
    bool v_found = false;
    size_t v = 0;
    while (!v_found && (v < xi.size())) {
        L = 0; U = xi[v].size() - 1;
        n->region(v, &L, &U);
        if (U>=L) v_found = true;
        v++;
    }
    return v_found;
}

// compute prob birth, goodbots save all bottom nodes that can be further splitted
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
    double pb;
    tree::npv bnv; //all bottom nodes
    t.getbots(bnv);
    for (size_t i=0; i!=bnv.size(); i++) {
        if (cansplit(bnv[i], xi))
            goodbots.push_back(bnv[i]);
    }
    if (goodbots.size()==0) {
        pb = 0.0;
    } else {
        if (t.treesize()==1) pb = 1.0;
        else pb = pi.pb;
    }
    return pb;
}

// find variables the node n can split on and store in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars)
{
    int L,U;
    for (size_t v=0; v!=xi.size(); v++) {
        L=0; U=xi[v].size()-1;
        n->region(v, &L, &U);
        if (U>=L) goodvars.push_back(v);
    }
}

// find variables the INTERNAL node n can split on and store in goodvars
void getinternalvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars)
{
    int L,U;
    for (size_t v=0; v!=xi.size(); v++) {
        L=0; U=xi[v].size()-1;
        getpertLU(n, v, xi, &L, &U);
        if (U>=L) goodvars.push_back(v);
    }
}

// number of avaible cut points for node n variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var)
{
    int L,U;
    
    getpertLU(n, var, xi, &L, &U);
    return std::max(0, U-L+1);
}

// get prob a node grows, 0 if no good vars, a/(1+d)^b else
double pgrow(tree::tree_p n, xinfo&xi, pinfo& pi)
{
    if (cansplit(n, xi)) {
        return pi.alpha/pow(1.0+n->depth(), pi.beta);
    } else {
        return 0.0;
    }
}

//calibart
// get the L,U values for a node in the tree GIVEN the tree structure below and above that node
void getLU(tree::tree_p pertnode, xinfo& xi, int* L, int* U)
{
    tree::tree_p l,r;
    
    *L=0;
    *U=xi[pertnode->getv()].size()-1;
    l=pertnode->getl();
    r=pertnode->getr();
    
    bool usel, user;
    usel = l->nuse(pertnode->getv());
    user = r->nuse(pertnode->getv());
    if (usel && user) {
        l->region_l(pertnode->getv(), L);
        r->region_u(pertnode->getv(), U);
    }
    else if (usel)
    {
        pertnode->region(pertnode->getv(), L, U);
        l->region_l(pertnode->getv(), L);
    }
    else
    {
        pertnode->region(pertnode->getv(), L, U);
        r->region_u(pertnode->getv(), U);
    }
}

/* Note: region is to search upwards based on current node, region_l is to search downwards knowing that the tree can be further splitted by rule v<=?, region_r is to search downwards and knowing that the tree can be further splitted by rule v>? */

// if the var is prescribed, given the nodes above and below it, what is the availale range for this node with this prescribes var
void getpertLU(tree::tree_p pertnode, size_t pertvar, xinfo& xi, int* L, int* U)
{
    *L=0;
    *U=xi[pertvar].size()-1;
    
    bool usel,user;
    usel=pertnode->l->nuse(pertvar);
    user=pertnode->r->nuse(pertvar);
    if (usel && user) {
        pertnode->l->region_l(pertvar, L);
        pertnode->r->region_u(pertvar, U);
    }
    else if (usel)
    {
        pertnode->region(pertvar, L, U);
        pertnode->l->region_l(pertvar, L);
    }
    else {
        pertnode->region(pertvar, L, U);
        pertnode->r->region_u(pertvar, U);
    }
}

//-------------------
// suff statistics for all bottom nodes
/* Note: Worker class is from RcppParallel*/
struct AllSuffWorker: public Worker
{
    tree& x;
    xinfo& xi;
    dinfo& di;
    size_t nb;
    std::map<tree::tree_cp, size_t> bnmap;
    double* weight;
    
    // internal state
    double n;
    double sy;
    double n0;
    
    std::vector<sinfo> sv_tmp;
    double *xx; //current x
    double y; //current y
    size_t ni; // index
    
    //constrcutor
    AllSuffWorker(tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp, size_t> bnmap, size_t nb, double* weight): x(x), xi(xi),di(di),nb(nb),bnmap(bnmap),weight(weight){
        n=0.0;
        sy=0.0;
        n0=0.0;
        sv_tmp.resize(nb);
    }
    
    AllSuffWorker(const AllSuffWorker& asw, Split):x(asw.x),xi(asw.xi),di(asw.di),nb(asw.nb),bnmap(asw.bnmap),weight(asw.weight){
        n=0.0;
        sy=0.0;
        n0=0.0;
        sv_tmp.resize(nb);
    }
    
    void operator()(std::size_t begin, std::size_t end){
        for (size_t i=begin; i<end; i++) {
            xx = di.x + i*di.p;
            y = di.y[i];
            
            ni = bnmap[x.bot(xx, xi)];
            
            sv_tmp[ni].n0 +=1;
            sv_tmp[ni].n += weight[i];
            sv_tmp[ni].sy += weight[i]*y;
        }
    }
    
    void join(const AllSuffWorker& asw){
        for (size_t i=0; i!=nb; i++) {
            sv_tmp[i].n0 += asw.sv_tmp[i].n0;
            sv_tmp[i].n += asw.sv_tmp[i].n;
            sv_tmp[i].sy += asw.sv_tmp[i].sy;
        }
    }
};

void allsuff(tree& x, xinfo& xi, dinfo& di, double* weight, tree::npv& bnv, std::vector<sinfo>& sv)
{
    tree::tree_cp tbn; // pointer to tree bottom mode for the current obs
    size_t ni; // index in vector of current bottom nodes
    double *xx;
    double y;
    
    bnv.clear();
    
    x.getbots(bnv);
    
    typedef tree::npv::size_type bvsz;
    bvsz nb = bnv.size();
    sv.resize(nb);
    
    std::map<tree::tree_cp, size_t> bnmap;
    for (bvsz i=0; i!=bnv.size(); i++) bnmap[bnv[i]]=i;
    AllSuffWorker asw(x,xi,di,bnmap,nb,weight);
    parallelReduce(0, di.N, asw);
    
    for (size_t i=0; i!=nb; i++) {
        sv[i].n0 += asw.sv_tmp[i].n0;
        sv[i].n += asw.sv_tmp[i].n;
        sv[i].sy += asw.sv_tmp[i].sy;
    }
}

// counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv)
{
    tree::tree_cp tbn;
    size_t ni;
    double *xx;
    double y;
    
    bnv.clear();
    x.getbots(bnv);
    
    typedef tree::npv::size_type bvsz;
    bvsz nb = bnv.size();
    
    std::vector<int> cts(bnv.size(),0);
    
    std::map<tree::tree_cp,size_t> bnmap;
    for (bvsz i=0; i!=bnv.size(); i++) bnmap[bnv[i]]=i;
    
    for (size_t i=0; i<di.N; i++) {
        xx = di.x + i*di.p; // vector for the ith subject
        y = di.y[i];
        
        tbn = x.bot(xx, xi); //find leaf for the ith subject;
        ni = bnmap[tbn]; // get its nid
        
        cts[ni] +=1; // count for this node, +1
    }
    return cts;
}
/*Note: cts is a vector of the number of obs within each leaf node of a tree*/

// update counts to reflect observation i?
void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, int sign)
{
    tree::tree_cp tbn;
    size_t ni;
    double *xx;
    double y;
    
    typedef tree::npv::size_type bvsz;
    bvsz nb = bnv.size();
    
    std::map<tree::tree_cp, size_t> bnmap;
    for (bvsz ii=1; ii!=bnv.size(); ii++) {
        bnmap[bnv[ii]]=ii;
    }
    
    xx = di.x + i*di.p;
    y = di.y[i];
    
    tbn = x.bot(xx, xi);
    ni = bnmap[tbn];
    
    cts[ni] += sign;
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp, size_t>& bnmap, int sign)
{
    tree::tree_cp tbn;
    size_t ni;
    double *xx;
    double y;
    // bnmap already given
    xx = di.x + i*di.p;
    y = di.y[i];
    
    tbn = x.bot(xx, xi);
    ni = bnmap[tbn];
    
    cts[ni] += sign;
    
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi, dinfo& di, std::map<tree::tree_cp, size_t>& bnmap, int sign, tree::tree_cp &tbn)
{
    // bnmap and tbn already exist;
    size_t ni;
    double *xx;
    double y;
    
    xx = di.x + i*di.p;
    y = di.y[i];
    
    tbn = x.bot(xx, xi);
    ni = bnmap[tbn];
    
    cts[ni] += sign;
}

// check min leaf size
bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di)
{
    bool good = true;
    tree::npv bnv;
    std::vector<int> cts;
    
    int m = 0;
    for (size_t tt=0; tt<t.size(); ++tt) {
        cts = counts(t[tt], xi, di, bnv);
        m = std::min(m, *std::min_element(cts.begin(), cts.end()));
        if (m<minct) {
            good = false;
            break;
        }
    }
    return good;
}

// ------------------------------
/*Note:skip all the MPI funcs for now*/
// ------------------------------

struct GetSuffBirthWorker: public Worker
{
    tree& x;
    tree::tree_cp nx;
    size_t v;
    size_t c;
    xinfo& xi;
    dinfo& di;
    double* phi;
    
    // internal state
    double l_n;
    double l_sy;
    double l_n0;
    double r_n;
    double r_n0;
    double r_sy;
    
    double *xx;
    double y;
    
    // Constructor
    
    GetSuffBirthWorker(tree &x,
                       tree::tree_cp nx,
                       size_t v,
                       size_t c,
                       xinfo& xi,
                       dinfo& di,
                       double* phi):x(x),nx(nx),v(v),c(c),xi(xi),di(di),phi(phi){
        
        l_n=0.0;
        l_sy=0.0;
        l_n0=0.0;
        
        r_n=0.0;
        r_sy=0.0;
        r_n0=0.0;
    }
    // splitting constructor
    GetSuffBirthWorker(const GetSuffBirthWorker& gsw, Split):x(gsw.x),nx(gsw.nx),v(gsw.v),c(gsw.c),xi(gsw.xi),di(gsw.di),phi(gsw.phi){
        
        l_n=0.0;
        l_sy=0.0;
        l_n0=0.0;
        
        r_n=0.0;
        r_sy=0.0;
        r_n0=0.0;
    }
    
    // An operator() which performs the work
    void operator()(std::size_t begin, std::size_t end){
        for (size_t i=begin; i<end; i++) {
            xx = di.x + i*di.p;
            if (nx == x.bot(xx, xi)) {
                y = di.y[i];
                if (xx[v] < xi[v][c]) {
                    l_n0 += 1;
                    l_n += phi[i]; //weights
                    l_sy += phi[i]*y;
                } else {
                    r_n0 += 1;
                    r_n += phi[i];
                    r_sy += phi[i]*y;
                }
            }
        }
    }
    
    void join(const GetSuffBirthWorker& gsw){
        l_n += gsw.l_n;
        l_sy += gsw.l_sy;
        l_n0 += gsw.l_n0;
        
        r_n += gsw.r_n;
        r_sy += gsw.r_sy;
        r_n0 += gsw.r_n0;
    }
};

// get sufficient stat for children (v,c) of node nx in tree x
void getsuffBirth(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr)
{
    GetSuffBirthWorker gsw(x,nx,v,c,xi,di,phi);
    
    parallelReduce(0,di.N,gsw);
    
    sl.n = gsw.l_n;
    sl.sy = gsw.l_sy;
    sl.n0 = gsw.l_n0;
    
    sr.n = gsw.r_n;
    sr.sy = gsw.l_sy;
    sr.n0 = gsw.r_n0;
}

struct GetSuffDeathWorker: public Worker
{
    tree& x;
    tree::tree_cp nl;
    tree::tree_cp nr;
    xinfo& xi;
    dinfo& di;
    double* phi;
    
    // internal state
    double l_n;
    double l_sy;
    double l_n0;
    
    double r_n;
    double r_sy;
    double r_n0;
    
    double *xx;
    double y;
    
    //constructor
    
    GetSuffDeathWorker(tree& x,
                       tree::tree_cp nl,
                       tree::tree_cp nr,
                       xinfo& xi,
                       dinfo& di,
                       double* phi):x(x),nl(nl),nr(nr),xi(xi),di(di),phi(phi){
        
        l_n = 0.0;
        l_sy = 0.0;
        l_n0 = 0.0;
        r_n = 0.0;
        r_sy = 0.0;
        r_n0 = 0.0;
    }
    
    //splitting constructor
    GetSuffDeathWorker(const GetSuffDeathWorker& gsw, Split):x(gsw.x),nl(gsw.nl),nr(gsw.nr),xi(gsw.xi),di(gsw.di),phi(gsw.phi){
        
        l_n = 0.0;
        l_sy = 0.0;
        l_n0 = 0.0;
        r_n = 0.0;
        r_sy = 0.0;
        r_n0 = 0.0;
    }
    
    void operator()(std::size_t begin, std::size_t end){
        for (size_t i=begin; i<end; i++) {
            xx = di.x + i*di.p;
            tree::tree_cp bn = x.bot(xx, xi);
            y = di.y[i];
            
            if (bn==nr) {
                r_n0 += 1;
                r_n += phi[i];
                r_sy += phi[i]*y;
            }
            if (bn == nl) {
                l_n0 += 1;
                l_n += phi[i];
                l_sy += phi[i]*y;
            }
        }
    }
    
    void join(const GetSuffDeathWorker& gsw){
        l_n += gsw.l_n;
        l_sy += gsw.l_sy;
        l_n0 += gsw.l_n0;
        
        r_n += gsw.r_n;
        r_sy += gsw.r_sy;
        r_n0 += gsw.r_n0;
    }
    
};

// get sufficient stat for pair of bottom children nl, nr, in tree x
void getsuffDeath(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr)
{
    GetSuffDeathWorker gsw(x,nl,nr,xi,di,phi);
    parallelReduce(0, di.N, gsw);
    
    sl.n = gsw.l_n;
    sl.sy = gsw.l_sy;
    sl.n0 = gsw.l_n0;
    
    sr.n = gsw.r_n;
    sr.sy = gsw.r_sy;
    sr.n0 = gsw.r_n0;
}

// ------------------------------
/*Note:skip all the MPI funcs for now*/
// ------------------------------
// log integrated likelihood, here mu_j of the node is integrated out.

double loglike(double n, double sy, double sigma, double tau)
{
    double d = 1/(tau*tau) + n; // n=\sum phi_i for heterogeneous?
    
    double out = -log(tau) - 0.5*log(d);
    out += 0.5*sy*sy/d;
    return out;
}
/*
struct FitWorker: public Worker
{
    tree& t;
    xinfo& xi;
    dinfo& di;
    
    // internal
    double *xx;
    tree::tree_cp bn;
    std::vector<double> &fv; //node means
    
    //constructor
    FitWorker(tree& t,
              xinfo& xi,
              dinfo& di,
              std::vector<double>& fv):t(t),xi(xi),di(di),fv(fv){ }
    
    // No split constructor?
    
    void operator()(std::size_t begin, std::size_t end){
        for (size_t i=begin; i<end; i++) {
            xx = di.x + i*di.p;
            bn = t.bot(xx, xi);
            fv[i] = bn->getmu();
        }
    }
    
    void join(const FitWorker& fir){ }
};

void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv)
{
    fv.resize(di.N);
    
    FitWorker fir(t,xi,di,fv);
    parallelReduce(0, di.N, fir);
}
*/
// without parallel work?
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
    double* xx;
    tree::tree_cp bn;
    for (size_t i=0; i<di.N; i++) {
        xx = di.x + i*di.p;
        bn = t.bot(xx, xi);
        fv[i] = bn->getmu();
    }
}

void predict(std::vector<tree>& tv, xinfo& xi, dinfo& di, double* fp)
{
    size_t ntemp = di.N;
    double* fptemp = new double[ntemp];
    
    for (size_t k=0; k<ntemp; k++) fp[k]=0.0;
    for (size_t j=0; j<tv.size(); j++) {
        fit(tv[j], xi, di, fptemp);
        for (size_t k=0; k<ntemp; k++) {
            fp[k] += fptemp[k];
        }
    }
    delete [] fptemp;
}

void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
    double *xx;
    tree::tree_cp bn;
    pv.resize(di.N);
    for (size_t i=0; i<di.N; i++) {
        xx = di.x + i*di.p;
        bn = t.bot(xx, xi);
        pv[i] = bn->nid();
    }
}

// draw all bottom nodes mu
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double* weight, RNG& gen)
{
    tree::npv bnv;
    std::vector<sinfo> sv;
    
    allsuff(t, xi, di, weight, bnv, sv);
    // bottom nodes stored in Expecting a single value: [extent=0].bnv
    // suff stat stores in sv, get from data
    // Do NOTE: sv[i].n is actually n/sigma^2 since weight = 1/sigma^2; and sv[i].y is thus \sum_y/sigma^2
    for (tree::npv::size_type i=0; i!=bnv.size(); i++) {
        double fcvar = 1.0/(1.0/(pi.tau*pi.tau)+sv[i].n);
        double fcmean = sv[i].sy*fcvar;
        bnv[i]->setmu(fcmean + gen.normal()*sqrt(fcvar));
        
        if (bnv[i]->getmu() != bnv[i]->getmu()) {
            for (int j=0; j<di.N; ++j) Rcout << *(di.x + j*di.p) << " "; // print covariate matrix
            Rcout << endl << "fcvar" << fcvar << "svi[n]" << sv[i].n << "i" << i;
            Rcout << endl << t;
            Rcpp::stop("drmu failed");
        }
    }
}

void prxi(xinfo& xi)
{
    Rcout << "xinfo:\n";
    for (size_t v=0; v!=xi.size(); v++) {
        Rcout << "v: " << v << endl;
        for (size_t j=0; j!=xi[v].size(); j++) {
            Rcout << "j,xi[v][j]: " << j << "," << xi[v][j] << endl;
        }
    }
    Rcout << "\n\n";
}

void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
    double xinc;
    
    std::vector<double> minx(p, INFINITY);
    std::vector<double> maxx(p, -INFINITY);
    double xx;
    for (size_t i=0; i<p; i++) {
        for (size_t j=0; j<n; j++) {
            xx = *(x+p*j+i);
            if (xx < minx[i]) minx[i]=xx;
            if (xx > maxx[i]) maxx[i]=xx;
        }
    }
    
    // make grid of nc cutpoints between min and max for each x
    xi.resize(p);
    for (size_t i=0; i<p; i++) {
        xinc = (maxx[i] - minx[i])/(nc+1.0);
        xi[i].resize(nc);
        for (size_t j=0; j<nc; j++) xi[i][j] = minx[i] + (j+1)*xinc;
    }
}

// get min/max for p predictors needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
    double xx;
    
    for (size_t i=0; i<p; i++) {
        for (size_t j=0; j<n; j++) {
            xx = *(x+p*j+i);
            if (xx < minx[i]) minx[i]=xx;
            if (xx > maxx[i]) maxx[i]=xx;
        }
    }
    
}

// make xinfo given min/max
void makexinfominmax(size_t p, xinfo&xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx)
{
    double xinc;
    xi.resize(p);
    for (size_t i=0; i<p; i++) {
        xinc = (maxx[i] - minx[i])/(nc+1.0);
        xi[i].resize(nc);
        for (size_t j=0; j<nc; j++) {
            xi[i][j] = minx[i] + (j+1)*xinc;
        }
    }
}

void updateLabels(int* labels, double* mixprop, double* locations, double* resid, double sigma, size_t nobs, size_t nmix, RNG& gen)
{
    // use random uniform generator to sample from multinomial
    double ptot, pcum, u;
    int count;
    for (size_t k=0; k<nobs; k++) {
        ptot = 0.0;
        count = 1;
        pcum = 0.0;
        for (size_t h=0; h<nmix; h++){
            ptot += mixprop[h]*(1/sqrt(2*PI*sigma*sigma))*exp(-0.5*(resid[k]-locations[h])*(resid[k]-locations[h])/(sigma*sigma));//R::dnorm(resid[k], locations[h], sigma, 0);
        }
        //u = ptot*R::runif(0.0, 1.0);
        u = ptot*gen.uniform();
        for (size_t h=0; h<nmix; h++) {
            pcum += mixprop[h]*(1/sqrt(2*PI*sigma*sigma))*exp(-0.5*(resid[k]-locations[h])*(resid[k]-locations[h])/(sigma*sigma));//R::dnorm(resid[k], locations[h], sigma, 0);
            if (u < pcum){
                break;
            }
            count += 1;
        }
        labels[k] = count; // Note: so the resulted count/labels will be 1 ~ nmix
    }
}

void updateMixprp(double* mixprop, double* mass, int* labels, int* mixcnts, size_t nmix, double psi1, double psi2, RNG& gen)
{
    int ncum;
    double shape1, shape2, gam_shape, gam_scale;
    double Vold, Vnew, log_mixprp, Vcum, tmp;
    
    ncum = std::accumulate(mixcnts+1, mixcnts+nmix, 0);
    shape1 = mixcnts[0] + 1.0;
    shape2 = mass[0] + ncum;
    Vold = gen.beta(shape1, shape2);
    //Vold = R::rbeta(shape1, shape2);
    log_mixprp = log(Vold);
    Vcum = log(1-Vold); // for updating mass
    mixprop[0] = Vold;
    for (size_t h=1; h<nmix-1; h++) {
        ncum = std::accumulate(mixcnts+h+1, mixcnts+nmix, 0);
        shape1 = mixcnts[h] + 1.0;
        shape2 = mass[0] + ncum;
        Vnew = gen.beta(shape1, shape2);
        //Vnew = R::rbeta(shape1, shape2);
        tmp = (Vnew/Vold)*(1-Vold);
        log_mixprp += log(tmp); // the log v_h
        
        Vcum += log1p(-Vnew); // for updating mass
        Vold = Vnew;
        mixprop[h] = exp(log_mixprp);
    }
    mixprop[nmix-1] = 1.0 - std::accumulate(mixprop, mixprop + nmix - 1, 0.0);
   
    gam_shape = psi1 + nmix - 1;
    gam_scale = 1/(psi2 - Vcum);

    //mass[0] = gen.gamma(psi1 + nmix - 1, 1/(psi2 - Vcum));
    mass[0] = R::rgamma(gam_shape, gam_scale);
   //Rprintf("Gamma draw': %f \n", Vcum);
}

void updateLocations(double* locations, double* mixprop, double* resid, int* labels, int* mixcnts, double sigma, double prior_sigsq, size_t nobs, size_t nmix, RNG& gen)
{
    double sigsq, compo_sum, compo_include, post_prec, post_mean, post_sd, muG;
    
    sigsq = sigma*sigma;
    for (size_t h=0; h<nmix; h++) {
        compo_sum = 0.0;
        for (size_t i=0; i<nobs; i++) {
            if (labels[i] == h+1) { // recall that the labels are 1,...,nmix
                compo_include = 1.0;
            } else {
                compo_include = 0.0;
            }
            //compo_include = (labels[i] == h+1) ? 1.0 : 0.0;
            compo_sum += compo_include*resid[i];
        }
        post_prec = 1/(prior_sigsq*mixcnts[h] + sigsq);
        post_mean = prior_sigsq*post_prec*compo_sum;
        post_sd = sqrt(sigsq*prior_sigsq*post_prec);
        locations[h] = gen.normal(post_mean, post_sd);
        
        //locations[h] = R::rnorm(post_mean, post_sd);
    }
    // normalization of the locations
    muG = 0.0;
    for (size_t h=0; h<nmix; h++) {
        muG += mixprop[h]*locations[h];
    }
    
    for (size_t h=0; h<nmix; h++) {
        locations[h] -= muG;
    }
}

void updateIndivLocations (double* indiv_locations, double* locations, int* labels, size_t nobs, size_t nmix)
{
    double loc = 0.0;
    for (size_t k=0; k<nobs; k++) {
        for (size_t h=0; h<nmix; h++) {
            if (labels[k] == h+1) {
                loc = locations[h];
                break;
            }
        }
        indiv_locations[k] = loc;
    }
}

double logPriT(tree::tree_p x, xinfo& xi, pinfo& pi)
{
    double p_grow = pgrow(x, xi, pi);
    double retval = 0.0;
    size_t v;
    std::vector<size_t> goodvars;
    int L,U;
    if (x->ntype() == 'b') {
        retval = log(1.0-p_grow);
    } else {
        retval = log(p_grow);
        getgoodvars(x, xi, goodvars);
        retval += log(1.0/(double)goodvars.size());
        v = x->getv();
        L=0;U=xi[v].size()-1;
        x->region(v, &L, &U);
        retval -= log((double)(U-L+1));
        retval += logPriT(x->getl(), xi, pi) + logPriT(x->getr(), xi, pi);
    }
    return retval;
}

bool CheckRule(tree::tree_p n, xinfo& xi, size_t var)
{
    int L,U;
    bool goodrule = false;
    L=0; U=xi[var].size()-1;
    n->region(var, &L, &U);
    if (!(n->getl())) {
        goodrule = true;
    } else {
        size_t v = n->getv();
        if (v == var) {
            size_t c = n->getc();
            if ((c>=L) && (c<=U)) {
                goodrule = (CheckRule(n->getl(), xi, var) && CheckRule(n->getl(), xi, var));
            }
        } else {
            goodrule = (CheckRule(n->getl(), xi, var) && CheckRule(n->getl(), xi, var));
        }
    }
    return goodrule;
}
