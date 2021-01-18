#ifndef MODELS_H
#define MODELS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

List eof(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
         arma::colvec m_0, arma::mat C_0, arma::mat scale_W,
         NumericVector params, CharacterVector proc_model,
         const int n_samples, const bool verbose);

List ide(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
         arma::mat mean_mu_kernel, arma::mat var_mu_kernel, arma::cube K,
         arma::cube scale_Sigma_kernel, arma::mat scale_W, NumericVector params, 
         const int n_samples, const bool verbose);

#endif
