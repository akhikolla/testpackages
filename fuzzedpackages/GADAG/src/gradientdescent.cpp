// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat gradientdescent (arma::mat P, int n, arma::mat XtX, float L, float lambda, int maxite, float tolobj){
  int p  = P.n_cols;
  arma::vec E = arma::zeros(maxite);
  arma::mat T0 = arma::zeros(p,p);
  
  arma::mat PtXtXP = P.t()*XtX*P;
  arma::mat I = arma::eye(p,p);
  arma::mat U = T0 + (2 * PtXtXP * (I-T0))/ (n*L);
  
  arma::mat Tnew = arma::sign(U) % max(arma::zeros(p,p),arma::abs(U)-lambda/L) - arma::trimatu(arma::sign(U) % max(arma::zeros(p,p),arma::abs(U)-lambda/L));
  
  int k = 1;
  float Tdiff = 1000000;
  while (k<(maxite+1) && Tdiff>tolobj){
    arma::mat Tk = Tnew;
    arma::mat Uk = Tk + (2 * PtXtXP * (I-Tk))/ (n*L);
    Tnew = arma::sign(Uk) % max(arma::zeros(p,p),arma::abs(Uk)-lambda/L) - arma::trimatu(arma::sign(Uk) % max(arma::zeros(p,p),arma::abs(Uk)-lambda/L));
    Tdiff = arma::accu((Tnew-Tk)%(Tnew-Tk));
    k = k+1;
  }
  return vectorise(Tnew,0).t();
}
