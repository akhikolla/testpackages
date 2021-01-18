#ifndef IDE_HELPERS_H
#define IDE_HELPERS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double kernelLikelihood(const arma::mat & G, const arma::mat & theta, 
                        const arma::mat & W);

double kernelLikelihoodDiscount(const arma::mat & G, const arma::mat & theta, 
                        const arma::cube & C, const double lambda);

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K);

arma::mat makeW(const int J, const double L);

arma::mat makeF(const arma::mat & locs, const arma::mat & w, const double L);

void makeB(arma::mat & B, const arma::mat & mu, const arma::cube & Sigma, 
           const arma::mat & locs, const arma::mat & w, const double L);

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K);

arma::mat proposeMu(arma::mat mu, arma::mat Sigma);

#endif
