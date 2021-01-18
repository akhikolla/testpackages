#ifndef SOFT_H
#define SOFT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double softc(const double &s, const double &tau);

void softmatrixc(arma::mat &S, const arma::mat &Tau);

int numzeros(arma::mat &X);

#endif
