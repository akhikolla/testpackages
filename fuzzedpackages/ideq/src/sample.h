#ifndef SAMPLE_H
#define SAMPLE_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void backwardSample(arma::mat & theta, const arma::mat & m, const arma::mat & a,
                    const arma::cube & C, const arma::mat & G, const arma::cube & R_inv);

void sampleAR(arma::mat & G, const arma::cube & W_inv, const arma::mat & theta,
              const arma::mat & Sigma_G_inv, const arma::mat & mu_G,
              const bool Discount = false, const double lambda = 0.0);

void sampleG(arma::mat & G, const arma::cube & W_inv, const arma::mat & theta,
             const arma::mat & Sigma_g_inv, const arma::mat & mu_g,
             const bool Discount = false, const double lambda = 0.0);

void sampleLambda(double & lambda_new, const double & alpha_lambda, const double & beta_lambda,
                  const arma::mat & G, const arma::cube & C, const arma::mat & theta);

void sampleSigma2(double & sigma2_new, const double & alpha_sigma2, const double & beta_sigma2,
                  const arma::mat & Y, const arma::mat & F, const arma::mat & theta);

void sampleW(arma::mat & W, const arma::mat & theta, const arma::mat & G,
             const arma::mat & scale_W, const int df_W);

#endif
