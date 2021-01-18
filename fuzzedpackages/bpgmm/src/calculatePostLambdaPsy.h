#ifndef __BPGMM_CALCULATEPOSTLAMBDAPSY__
#define __BPGMM_CALCULATEPOSTLAMBDAPSY__

#include <RcppArmadillo.h>

Rcpp::List Calculate_PostLambdaPsy(Rcpp::S4 hparam, Rcpp::List CxyList, Rcpp::S4 thetaYList, arma::vec constraint);


#endif
