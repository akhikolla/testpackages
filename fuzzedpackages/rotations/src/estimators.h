#ifndef ESTIMATORS_H
#define ESTIMATORS_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
int checkQ4(const Rcpp::NumericMatrix &Q);

// [[Rcpp::export]]
int checkSO3(const arma::mat &Rs);

// [[Rcpp::export]]
arma::mat expskewC(const arma::mat &M);

// [[Rcpp::export]]
arma::mat expskewCMulti(const arma::mat &M);

// [[Rcpp::export]]
arma::mat logSO3C(const arma::mat &R);

// [[Rcpp::export]]
arma::mat logSO3CMulti(const arma::mat &R);

// [[Rcpp::export]]
arma::mat projectSO3C(const arma::mat &M);

// [[Rcpp::export]]
arma::mat meanSO3C(const arma::mat &Rs);

// [[Rcpp::export]]
arma::rowvec meanQ4C(const arma::mat &Q);

// [[Rcpp::export]]
arma::mat medianSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps);

//[[Rcpp::export]]
arma::mat HartmedianSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps);

// [[Rcpp::export]]
arma::mat gmeanSO3C(const arma::mat &Rs, unsigned int maxIterations, double maxEps);

#endif /* ESTIMATORS_H */
