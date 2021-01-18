#include <RcppArmadillo.h>

RcppExport SEXP getInitial( SEXP residual_, SEXP phi_, SEXP S_){
  Rcpp::NumericVector residual(residual_);
  Rcpp::NumericMatrix phi(phi_); //Create Rcpp matrix from SEXP
  Rcpp::NumericVector standres;
  arma::mat PHI(phi.begin(), phi.nrow(), phi.ncol(), false); // reuses memory and avoids extra copy
  arma::mat ttX;
  int S = Rcpp::as<int>(S_);
  double sxy = 0, sxx = 0;
  double prop, sig2, sig2g;

  double sig2_sig2g = Rcpp::sum(Rcpp::pow(residual,2)) / (S - 1);
  standres = residual / sqrt(sig2_sig2g);

  ttX = PHI - arma::eye<arma::mat>(S, S);
  ttX = ttX - arma::mean(arma::mean(ttX));

  for(int i = 0; i < S - 1; ++i){
    for(int j = 1; j < S; ++j){
      sxy = sxy + standres[i] * standres[j] * ttX(i, j);
      sxx = sxx + pow(ttX(i, j), 2);
    }
  }

  prop = sxy / sxx;

  if(prop < 0)
    prop = 0;

  sig2g = prop *sig2_sig2g;
  sig2 = sig2_sig2g - sig2g;

  return Rcpp::List::create(Rcpp::Named("sig2") = sig2, Rcpp::Named("sig2g") = sig2g);
}

