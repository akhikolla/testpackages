#include <RcppArmadillo.h>

//Function to get XTheta matrix-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetXTheta(arma::vec const& Theta, arma::uvec const& XThetaInd, arma::vec const& TimeVec, arma::vec const& OneNu, arma::vec const& OneN, double tNu, int N, int M) {
  arma::vec ThetaLong = arma::kron(OneNu, Theta);
  arma::mat Stacks(N, 2, arma::fill::ones), XTheta(N, 2 * M, arma::fill::zeros);
  Stacks.col(1) = (TimeVec - ThetaLong) % (1 * (ThetaLong <= TimeVec));
  XTheta(XThetaInd) = Stacks;
  return XTheta;
}



//Function to get XTheta matrix at a particular location-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetXThetaLoc(double ThetaLoc, arma::vec const& Time, arma::vec const& OneNu, int Nu) {
  arma::vec ThetaLong = ThetaLoc * OneNu;
  arma::mat XThetaLoc(Nu, 2, arma::fill::ones);
  XThetaLoc.col(1) = (Time - ThetaLong) % (1 * (ThetaLong <= Time));
  return XThetaLoc;
}



//Create phi from random effects-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec CreatePhi(arma::vec const& Beta, arma::vec const& Lambda, arma::vec const& Eta, int M) {
  arma::mat PhiMatrix(5, M);
  arma::uvec b(2); b(0) = 0, b(1) = 1;
  arma::uvec l(2); l(0) = 2, l(1) = 3;
  arma::uvec e(1); e(0) = 4;
  PhiMatrix.rows(b) = arma::reshape(Beta, 2, M);
  PhiMatrix.rows(l) = arma::reshape(Lambda, 2, M);
  PhiMatrix.rows(e) = arma::trans(Eta);
  return arma::vectorise(PhiMatrix);
}



//Matrix inverse using cholesky decomoposition for covariances-------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CholInv(arma::mat const& Cov) {
  return arma::inv_sympd(Cov);
}



//Matrix inverse of for 3x3 matrix-----------------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Inv3(arma::mat const& A) {
  arma::mat result = arma::mat(3, 3);
  double determinant = A(0, 0) * ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) - A(0, 1) * ( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) + A(0, 2) * ( A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0) );
  double invdet = 1 / determinant;
  result(0, 0) =  ( A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2) ) * invdet;
  result(1, 0) = -( A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1) ) * invdet;
  result(2, 0) =  ( A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1) ) * invdet;
  result(0, 1) = -( A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0) ) * invdet;
  result(1, 1) =  ( A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) ) * invdet;
  result(2, 1) = -( A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2) ) * invdet;
  result(0, 2) =  ( A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1) ) * invdet;
  result(1, 2) = -( A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1) ) * invdet;
  result(2, 2) =  ( A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1) ) * invdet;
  return result;
}



//Matrix inverse of for 2x2 matrix-----------------------------------------------------------------------
arma::mat Inv2(arma::mat const& A) {
  arma::mat result = arma::mat(2, 2);
  double determinant = ( A(0, 0) * A(1, 1) ) - ( A(0, 1) * A(0, 1) );
  double invdet = 1 / determinant;
  result(0, 0) =  A(1, 1) * invdet;
  result(0, 1) = -A(1, 0) * invdet;
  result(1, 0) = -A(0, 1) * invdet;
  result(1, 1) =  A(0, 0) * invdet;
  return result;
}



//Function for making an upper diagonal matrix symmetric-------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat makeSymm(arma::mat const& A) {
  return arma::symmatu(A);
}



//Function that checks numerical equality of two objects against a tolerance-----------------------------
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}
