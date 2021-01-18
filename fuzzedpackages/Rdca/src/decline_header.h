
#ifndef DECLINE_HEADER_H
#define DECLINE_HEADER_H

#include "RcppArmadillo.h"

arma::mat exponential(Rcpp::List lst, arma::vec time);
arma::mat harmonic(Rcpp::List lst, arma::vec time);
arma::mat hyperbolic(Rcpp::List lst, arma::vec time);
arma::mat modified_hyperbolic(Rcpp::List lst, arma::vec time);

#endif
