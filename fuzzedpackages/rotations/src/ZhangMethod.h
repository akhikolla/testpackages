#ifndef ZHANGMETHOD_H
#define ZHANGMETHOD_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
Rcpp::NumericVector RdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2);

// [[Rcpp::export]]
arma::rowvec rdistSO3C(const arma::mat &Rs, const arma::mat &R2);

// [[Rcpp::export]]
Rcpp::NumericVector EdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2);

// [[Rcpp::export]]
double oneRdistC(const Rcpp::NumericMatrix &Q1, const Rcpp::NumericVector &Q2);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsC(const Rcpp::NumericMatrix &Qs, const Rcpp::NumericVector &Qhat);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsCMedian(
    const Rcpp::NumericMatrix &Qs,
    const Rcpp::NumericVector &Qhat
);

// [[Rcpp::export]]
Rcpp::NumericVector zhangQ4(const Rcpp::NumericMatrix &Q, unsigned int m);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsCSO3(const arma::mat &Rs, const arma::mat &Rhat);

// [[Rcpp::export]]
Rcpp::NumericVector zhangMedianC(const arma::mat &Rs, unsigned int m);

#endif /* ZHANGMETHOD_H */
