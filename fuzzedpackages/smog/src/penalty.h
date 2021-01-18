#ifndef PENALTY_H
#define PENALTY_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// proximate on L1 penalty
double proxL1(const double &x, const double &lambda);

// proximate on L2 penalty
arma::vec proxL2(const arma::vec &s, const double &lambda);

// proximate operators for different hierarchies
arma::vec prox(const arma::vec &s, const arma::vec &lambda, 
               const int &hierarchy, const arma::uvec &d);

// penalty function for different hierarchies
double penalty(const arma::vec &x, const arma::vec &lambda, 
               const int &hierarchy, const arma::uvec &d);

Rcpp::List glog(const arma::mat & y, 
                const arma::mat & x, 
                const arma::uvec & g, 
                const arma::uvec & v,
                const arma::vec & lambda,
                const int & hierarchy,
                const std::string & family,
                const double & rho,
                const bool & scale,
                const double & eabs,
                const double & erel,
                const double & LL,
                const double & eta,
                const int & maxitr);

#endif
