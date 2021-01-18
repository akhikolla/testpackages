// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// computes the spatial covariance matrix for a given kernel matrix using sparse matrix
//
// [[Rcpp::export]]
arma::mat sp_lcov_sparse(const arma::mat & x,
                         const arma::mat & k) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::mat L(p,p);
  L.zeros();
  
  arma::sp_mat k_sp = arma::sp_mat(k);
  
  arma::sp_mat::iterator it     = k_sp.begin();
  arma::sp_mat::iterator it_end = k_sp.end();
  
  for(; it != it_end; ++it)
  {
    L = L + (*it) * (x.row(it.row()).t()*x.row(it.col()));
  }  
  
  L = L/n;
  
  return L;
}

// computes the local difference matrix for a given kernel matrix using sparse matrix
//
// [[Rcpp::export]]
arma::mat sp_ldiff_sparse(const arma::mat & x,
                          const arma::mat & k) {
  const int n = x.n_rows;
  const int p = x.n_cols;
  arma::rowvec diff(p);
  diff.zeros();
  arma::mat L(p,p);
  L.zeros();
  
  arma::sp_mat k_sp = arma::sp_mat(k);
  
  arma::sp_mat::iterator it     = k_sp.begin();
  arma::sp_mat::iterator it_end = k_sp.end();
  
  for(; it != it_end; ++it) {
    diff = x.row(it.row()) - x.row(it.col());
    L += (*it) * (diff.t() * diff);
  }  
  
  L = L / n;
  
  return L;
}


