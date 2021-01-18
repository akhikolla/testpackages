#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

arma::mat MultInv(
  const arma::mat& MatrixToInverse
) {
  mat ResultMatrix;

  // Try inv_sympd
  bool inv_sympd_success = inv_sympd(ResultMatrix, MatrixToInverse);

  // If unsuccessful, try inv
  if (inv_sympd_success == FALSE) {
    ResultMatrix = inv(MatrixToInverse);
  }

  // Either one success or one error thrown
  return ResultMatrix;
}
