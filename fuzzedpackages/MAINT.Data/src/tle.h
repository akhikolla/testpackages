#ifndef TLEH
#define TLEH

#include <vector>
//#include <Rcpp.h>
//#include <RcppEigen.h>
#include "RcppArmadillo.h"
#include <numeric>

using namespace std;
using namespace Rcpp;

const double Inf = DBL_MAX;
void parcolmeans(const NumericMatrix& X,const vector<int>& Set,arma::vec& res);

class Estimate {
  public: 
//    Estimate(int p) : p_(p) { muE_.resize(p); SigmaE_.resize(p,p); };
    Estimate(int p) : p_(p) { muE_.resize(p); SigmaE_.resize(p,p); };
    int p(void) { return p_; }
    arma::vec& muE(void) { return muE_; }
    arma::mat& SigmaE(void) { return SigmaE_; }
    double logLik(void) { return logLik_; }
    void setlogLik(double LogL) { logLik_ = LogL; } 
  private:
    int p_;
    arma::vec muE_;
    arma::mat SigmaE_;
    double logLik_;		
};

RcppExport SEXP Cfasttle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP Poolm_s, SEXP m_s, SEXP kdblstar_s, SEXP k_s, SEXP nrep_s,
  SEXP Cnf_s, SEXP c0_s, SEXP maxrefstps_s, SEXP limlnk2_s, SEXP ClctSt_s);

RcppExport SEXP Cfulltle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP k_s, SEXP Cnf_s, SEXP limlnk2_s, SEXP c0_s);



#endif
