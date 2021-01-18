// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// computes idw interpolation
//
// [[Rcpp::export]]
arma::mat idw(const arma::mat& coords_pred,
              const arma::mat& coords_vals,
              const arma::mat& vals,
              const int& p) {
  
  const int n = coords_pred.n_rows;
  const int m = coords_vals.n_rows;
  
  double distance;

  arma::mat w_mat = arma::zeros<arma::mat>(n, m);
  
  for(int i=0; i < n; ++i) {
    for(int j=0; j < m; ++j) {
      distance = norm(coords_pred.row(i) - coords_vals.row(j));
      if(distance == 0){
        w_mat.row(i) = arma::zeros<arma::mat>(1, m); 
        w_mat(i,j) = 1;
        break;
      } else {
        w_mat(i, j) = 1.0 / pow(distance, p);
      }
    }
  }
  
  arma::mat w_sum = arma::zeros<arma::mat>(n, n);
  w_sum.diag() = pow(sum(w_mat, 1), -1);
  
  w_mat = w_sum * w_mat * vals;
  
  return w_mat;
}
