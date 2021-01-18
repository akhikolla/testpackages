#ifndef PARAMS_H
#define PARAMS_H

#include <vector>
#include <RcppArmadillo.h>  


/* DYNAMIC CREATION OF LOCAL VECTORS AND MATRICES */
typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<IntVec> IntMatrix;
typedef std::vector<DoubleVec> DoubleMatrix;

using namespace std;
using namespace Rcpp;
using namespace arma;

class Params { // maybe include beta someday?
public:
  NumericVector w;
  int nn; // number of locations
  int nn2; // number of observations, including multiple obs at same loc
  int kk;
  double p0; // or alpha
  int nreg; // number of observations expected per block for regression model
  double d;
  mat sigma_jitter;
  DoubleVec priors;
  bool reg;
  bool graph;
  
  // constructor
  Params(SEXP pw, int pnn, int pnn2, SEXP pp0, bool preg, bool pgraph,
         double pba = 0, int pnreg = 0, int nDim = 1) {
    reg = preg;
    graph = pgraph;
    w = NumericVector(pw); // w2
    nn = pnn;
    nn2 = pnn2;
    p0 = NUMERIC_DATA(pp0)[0];
    if (reg) { // regression case
      kk = w.size()-1;
      sigma_jitter = ones(kk,kk)*0.01;
      d = pba;
      nreg = pnreg;
      
      double tmp;
      for (int j = 1; j < nn - 2; j++) { // maybe it can be bigger?
        tmp = Rf_pbeta(p0, (double) j, (double) nn - j + 1, 1, 1)
        + Rf_lbeta((double) j, (double) nn - j + 1);
        priors.push_back(tmp);
      }
    } else {
      // multivariate case
      kk = nDim;
      // a temporary bandaid for the large sample size cases
      double tmp;
      for (int j = 1; j < nn - 2; j++) { // maybe it can be bigger?
        tmp = exp(Rf_lbeta((double) j+1, (double) nn - j)) * 
               Rf_pbeta(p0, (double) j+1, (double) nn - j, 1, 0)/
              (exp(Rf_lbeta((double) j, (double) nn - j + 1)) * 
                Rf_pbeta(p0, (double) j, (double) nn - j + 1, 1, 0));
        priors.push_back(tmp);
      }
    } 
    

  }
  
  //methods
  void print() {
    Rprintf("Params-- nn: %d kk %d p0:%0.2f d:%0.2f nreg:%d\n", 
            nn, kk, p0, d, nreg);
  }
};


class GraphParams: public Params {
public:
  double p1; // probability of correct APP
  int freqAPP;
  int boundaryType; // type = 1 is node, type = 2 is edge
  
  // some MCMC parameters (probably shouldn't go here, but let's stick with it
  // for convenience)
  int burnin;
  int mcmc;
  bool doneBurnin;
  GraphParams(int nrow_y, SEXP pw, SEXP pa, int numLocs, SEXP type, 
              SEXP pburnin, SEXP pmcmc, SEXP pp1, SEXP pfreqAPP,
              SEXP pba, SEXP pnreg): Params(pw, numLocs, numLocs, pa, true, true,
              NUMERIC_DATA(pba)[0], INTEGER_DATA(pnreg)[0]) { // for regression
    nn2 = nrow_y;
    kk = w.size()-1;
    sigma_jitter = ones(kk,kk)*0.01;
    
    boundaryType = INTEGER_DATA(type)[0];
    burnin = INTEGER_DATA(pburnin)[0];
    mcmc = INTEGER_DATA(pmcmc)[0];
    p1 = NUMERIC_DATA(pp1)[0];
    freqAPP = INTEGER_DATA(pfreqAPP)[0];
    
    doneBurnin = false;
  }
  
  // yup a lot of variables are undeclared... but that's okay.
  GraphParams(SEXP pw, SEXP pa, int numLocs, int numNodes, SEXP type, SEXP pburnin,
              SEXP pmcmc, SEXP pp1, int pmm) : 
    Params(pw, numLocs, numNodes, pa, false, true, 0, 0, pmm)
  {
    boundaryType = INTEGER_DATA(type)[0];
    burnin = INTEGER_DATA(pburnin)[0];
    mcmc = INTEGER_DATA(pmcmc)[0];
    p1 = NUMERIC_DATA(pp1)[0];
  }
  void print() {
    if (!reg) {
      // not regression
      Rprintf("multivariate: alpha=%0.2f, d=%0.2f w:%0.2f kk:%d\n", p0, d, w[0], kk);
    } else {
      // regression since multiple predictors
      Rprintf("regression: alpha=%0.2f, d=%0.2f w.size:%d kk:%d\n", p0, d, w.size(), kk);
    }
  }
  
};


#endif
