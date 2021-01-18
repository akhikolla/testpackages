#include <RcppArmadillo.h>

//Calculate inverse root of upper cholesky for log density of normal-----------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye) {
  return arma::solve(arma::trimatu(arma::chol(Cov)), Eye);
}



//Create function for creating a cube of covariance matrices at each visit------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube JointCovarianceCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu) {
  arma::cube JointCovariances(M, M, Nu);
  arma::mat Q(M, M), WStar(M, M), WAlpha(M, M), DWAlpha(M, M, arma::fill::zeros);
  for (int i = 0; i < Nu; i++) {
    WAlpha = WAlphas.slice(i);
    DWAlpha.diag() = arma::sum(WAlpha, 1);
    WStar = DWAlpha - WAlpha;
    Q = Rho * WStar + (1 - Rho) * EyeM;
    JointCovariances.slice(i) = Tau2(i) * arma::inv_sympd(Q);
  }
  return JointCovariances;
}



//Function for computing tau2_t Q(alpha_t)^-1 for one time point---------------------------------------------------
arma::mat JointCovarianceMatrix(arma::mat const& WAlpha, double tau2, arma::mat const& EyeM, double Rho, int M) {
  arma::mat Q(M, M), WStar(M, M), DWAlpha(M, M, arma::fill::zeros);
  DWAlpha.diag() = arma::sum(WAlpha, 1);
  WStar = DWAlpha - WAlpha;
  Q = Rho * WStar + (1 - Rho) * EyeM;
  return tau2 * arma::inv_sympd(Q);
}



//Create function for creating a cube of precision matrices at each visit------------------------------------------
arma::cube JointPrecisionCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu) {
  arma::cube JointPrecisions(M, M, Nu);
  arma::mat Q(M, M), WStar(M, M), WAlpha(M, M), DWAlpha(M, M, arma::fill::zeros);
  for (int i = 0; i < Nu; i++) {
    WAlpha = WAlphas.slice(i);
    DWAlpha.diag() = arma::sum(WAlpha, 1);
    WStar = DWAlpha - WAlpha;
    Q = Rho * WStar + (1 - Rho) * EyeM;
    JointPrecisions.slice(i) = Q / Tau2(i);
  }
  return JointPrecisions;
}



//Function to create array of Rooti likelihood matrices------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube RootiLikelihoodCube(arma::cube const& JointCovariances, arma::mat const& EyeM, int M, int Nu) {
  arma::cube RootiLikelihoods(M, M, Nu);
  for (int i = 0; i < Nu; i++) {
    RootiLikelihoods.slice(i) = GetRooti(JointCovariances.slice(i), EyeM);
  }
  return RootiLikelihoods;
}



//Function to calculate temporal correlation structure-------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SIGMA(double Phi, int TempCorInd, arma::mat const& TimeDist, int Nu) {

  //Create ouput object
  arma::mat out(Nu, Nu);

  //exponential
  if (TempCorInd == 0) out = arma::exp(-Phi * TimeDist);

  //ar(1) continuous
  if (TempCorInd == 1) {
    out = arma::eye(Nu, Nu);
    for (int j = 0; j < Nu; j++) {
      for (int k = 0; k < j; k++) {
        out(j, k) = std::pow(Phi, TimeDist(j, k));
      }
    }
    out = arma::symmatl(out);
  }
  return out;
}



//Create function for creating a cube of W(alpha_t) at each visit--------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube WAlphaCube(arma::vec const& Alpha, arma::colvec const& Z, arma::mat const& W, int M, int Nu, int WeightsInd) {
  arma::mat UpperW = arma::trimatu(W);
  arma::umat AdjacentEdgesBoolean = find(UpperW == 1);
  arma::cube WAlphas(M, M, Nu);
  arma::mat WAdj(M, M, arma::fill::zeros);
  int Zlength = Z.n_rows;
  arma::colvec Temp(Zlength);
  double alpha;
  for (int i = 0; i < Nu; i++) {
    alpha = Alpha(i);
    Temp = arma::exp(-Z * alpha);
    if (WeightsInd == 0) WAdj(AdjacentEdgesBoolean) = Temp; //continuous weights
    if (WeightsInd == 1) { //binary weights from Lee and Mitchell 2011
      Temp.elem(find(Temp >= 0.5)).ones();
      Temp.elem(find(Temp < 0.5)).zeros();
      WAdj(AdjacentEdgesBoolean) = Temp;
    }
    WAdj = arma::symmatu(WAdj);
    WAlphas.slice(i) = WAdj;
  }
  return WAlphas;
}



//Function to calculate W(alpha_t)---------------------------------------------------------------------------------
arma::mat WAlphaMatrix(double alpha, arma::colvec const& Z, arma::umat const& AdjacentEdgesBoolean, arma::mat const& W, int M, int WeightsInd) {
  arma::mat WAdj(M, M, arma::fill::zeros);
  arma::colvec Temp = arma::exp(-Z * alpha);
  if (WeightsInd == 0) WAdj(AdjacentEdgesBoolean) = Temp; //continuous weights
  if (WeightsInd == 1) { //binary weights from Lee and Mitchell 2011
    Temp.elem(find(Temp >= 0.5)).ones();
    Temp.elem(find(Temp < 0.5)).zeros();
    WAdj(AdjacentEdgesBoolean) = Temp;
  }
  WAdj = arma::symmatu(WAdj);
  return WAdj;
}
