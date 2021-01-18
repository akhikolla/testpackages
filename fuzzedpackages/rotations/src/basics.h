#ifndef BASICS_H
#define BASICS_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat eskewC(const arma::rowvec &U);

// [[Rcpp::export]]
arma::mat SO3defaultC(const arma::mat &U, const arma::vec &theta);

// [[Rcpp::export]]
arma::mat Q4defaultC(const arma::mat &U, const arma::vec &theta);

// [[Rcpp::export]]
arma::mat pMatC(const arma::mat &p);

// [[Rcpp::export]]
arma::mat genrC(const arma::vec &r, const arma::mat &S, int SO3, const arma::mat &u);

#endif /* BASICS_H */
