#ifndef ROBUSTCPP_H
#define ROBUSTCPP_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::rowvec HnCpp(const arma::mat &Qs);

arma::rowvec RdistCArma(const arma::mat &Q1, const arma::rowvec &Q2);

// [[Rcpp::export]]
arma::rowvec HnCppIntrinsic(const arma::mat &Qs);

// [[Rcpp::export]]
arma::rowvec HnCppBloc(const arma::mat &Qs, const arma::mat &Cs);

#endif /* ROBUSTCPP_H */
