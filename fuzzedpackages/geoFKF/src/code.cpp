#include <RcppArmadillo.h>

//' @name logLikMultiNorm
//' @title Log likelihood function for multivariate normal with spatial dependency.
//'
//' @param mCoef coefficient matrix. Each column is the coefficient from a curve;
//' @param mDist distance matris;
//' @param s2 variance from the covariance model;
//' @param phi variance from the covariance model;
//' @param rho variance from the covariance model;
//' @export

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::colvec logLikMultiNorm(const arma::mat& mCoef, const arma::mat& mDist, double s2, double phi, double p) {
  int nLoc = mCoef.n_cols; // number of points
  int nCoef = mCoef.n_rows; // size of vector of coefficients
  arma::colvec x = vectorise(mCoef, 0) - kron(arma::ones(nLoc, 1), mean(mCoef, 1));
  arma::mat mAR = arma::ones(nCoef, nCoef);
  arma::mat mCov = s2 * exp(- mDist / phi);
  arma::colvec output = arma::zeros(1);

  // matrix of AR part of covariance matrix
  for(int i = 0; i < nCoef - 1; i++) {
    for(int j = i + 1; j < nCoef; j++) {
      mAR(i, j) = pow(p, j - i);
      mAR(j, i) = pow(p, j - i);
    }
  }

  output = 0.5 * trans(x) * kron(inv(mCov), inv(mAR)) * x + 0.5 * nLoc * log(det(mAR)) + 0.5 * nCoef * log(det(mCov));

  return output;
}
