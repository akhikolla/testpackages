#ifndef __built_in_distributions_h__
#define __built_in_distributions_h__

#include <RcppArmadillo.h>

double sWeibull(double t, const Rcpp::List distPars);
double sBurr(double t, const Rcpp::List distPars);
double sgamma(double t, const Rcpp::List distPars);
double sgengamma(double t, const Rcpp::List distPars);
double surv (double t, const Rcpp::List distPars, const std::string dist);
arma::vec getextrapolPars(const Rcpp::List distPars, const std::string dist);

#endif 
