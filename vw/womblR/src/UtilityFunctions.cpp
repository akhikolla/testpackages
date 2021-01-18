#include <RcppArmadillo.h>

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
