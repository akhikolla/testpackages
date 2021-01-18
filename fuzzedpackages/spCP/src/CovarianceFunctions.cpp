#include <RcppArmadillo.h>

//Calculate inverse root of upper cholesky for log density of normal-----------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye) {
  return arma::solve(arma::trimatu(arma::chol(Cov)), Eye);
}



//Function for computing Q(alpha)^-1 for one time point---------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat QInvFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M) {
  arma::mat Q(M, M), WStar(M, M), DWAlpha(M, M, arma::fill::zeros);
  DWAlpha.diag() = arma::sum(WAlpha, 1);
  WStar = DWAlpha - WAlpha;
  Q = Rho * WStar + (1 - Rho) * EyeM;
  return arma::inv_sympd(Q);
}



//Function for computing Q(alpha) for one time point---------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat QFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M) {
  arma::mat Q(M, M), WStar(M, M), DWAlpha(M, M, arma::fill::zeros);
  DWAlpha.diag() = arma::sum(WAlpha, 1);
  WStar = DWAlpha - WAlpha;
  Q = Rho * WStar + (1 - Rho) * EyeM;
  return Q;
}


//Function to calculate W(alpha)---------------------------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat WAlphaFnc(double Alpha, arma::colvec const& DMLong, arma::umat const& AdjacentEdgesBoolean, arma::Mat<int> const& W, int M, int WeightsInd) {
  arma::mat WAdj(M, M, arma::fill::zeros);
  arma::colvec Temp = arma::exp(-DMLong * Alpha);
  if (WeightsInd == 0) WAdj(AdjacentEdgesBoolean) = Temp; //continuous weights
  if (WeightsInd == 1) { //binary weights from Lee and Mitchell 2011
    Temp.elem(find(Temp >= 0.5)).ones();
    Temp.elem(find(Temp < 0.5)).zeros();
    WAdj(AdjacentEdgesBoolean) = Temp;
  }
  WAdj = arma::symmatu(WAdj);
  return WAdj;
}
