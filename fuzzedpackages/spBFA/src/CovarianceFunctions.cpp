#include <RcppArmadillo.h>

//Calculate inverse root of upper cholesky for log density of normal-----------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye) {
  return arma::solve(arma::trimatu(arma::chol(Cov)), Eye);
}



//Function to calculate temporal correlation structure-------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat H(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu) {

  //Create ouput object
  arma::mat out(Nu, Nu);

  //exponential
  if (TempCorInd == 0) out = arma::exp(-Psi * TimeDist);

  //ar(1) continuous
  if (TempCorInd == 1) {
    out = arma::eye(Nu, Nu);
    for (int j = 0; j < Nu; j++) {
      for (int k = 0; k < j; k++) {
        out(j, k) = std::pow(Psi, TimeDist(j, k));
      }
    }
    out = arma::symmatl(out);
  }
  return out;
}



//Function to calculate spatial correlation structure-------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SpEXP(double Rho, arma::mat const& SpDist, int M) {
  arma::mat out(M, M);
  out = arma::exp(-Rho * SpDist);
  return out;
}
