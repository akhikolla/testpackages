#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

void MuVar(
  long double* return_var,
  const arma::mat& Predictors,
  const arma::mat& muCov,
  const arma::mat& muVarCov
) {
  const int N_rows = Predictors.n_rows;

  // Calculating the column sums
  const rowvec colsums = sum(Predictors);

  // Storing the mu-estimation and the mu-estimation-variance-estimation in
  // the provided return variable.
  *return_var = as_scalar(colsums * muCov) / N_rows;
  *(return_var + 1) = as_scalar(colsums * muVarCov * colsums.t())
    / (N_rows * N_rows);
}
