#ifndef SIGMA_H
#define SIGMA_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat RIDGEc(const arma::mat &S, double lam);

Rcpp::List ADMMc(const arma::mat &S, const arma::mat &initOmega, const arma::mat &initZ, const arma::mat &initY, const double lam, const double alpha = 1, bool diagonal = false, double rho = 2, const double mu = 10, const double tau_inc = 2, const double tau_dec = 2, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int maxit = 1e4);

#endif
