// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// computes the kernel matrix for a ball kernel
//
// [[Rcpp::export]]
arma::mat k_mat_ball(const arma::mat & coords,
                     const double & h) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.ones();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      if(distance > h){
        k(i,j) = k(j,i) = 0;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a ring kernel
//
// [[Rcpp::export]]
arma::mat k_mat_ring(const arma::mat & coords,
                     const double & h1,
                     const double & h2) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.zeros();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      if(h1 < distance && distance <= h2){
        k(i,j) = k(j,i) = 1;
      }
    }
  }
  
  return k;
}

// computes the kernel matrix for a gauss kernel
//
// [[Rcpp::export]]
arma::mat k_mat_exp(const arma::mat & coords,
                    const double & h) {
  const int n = coords.n_rows;
  double distance;
  arma::mat k(n, n);
  k.ones();
  
  for(int i=0; i < n; ++i) {
    for(int j=i+1; j < n; ++j) {
      distance = norm(coords.row(i) - coords.row(j));
      k(i,j) = k(j,i) = exp(- 0.5 * pow(distance, 2) / h);
    }
  }
  
  return k;
}


