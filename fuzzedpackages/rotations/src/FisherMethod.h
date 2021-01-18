#ifndef FISHERMETHOD_H
#define FISHERMETHOD_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
double fisherAxisC(const arma::mat &Qs, const arma::rowvec &Qhat);

// [[Rcpp::export]]
double fisherAxisCSymmetric(const arma::mat &Qs, const arma::rowvec &Qhat);

// [[Rcpp::export]]
arma::vec fisherBootC(const arma::mat &Qs, unsigned int m, bool symm);

#endif /* FISHERMETHOD_H */
