#include <RcppArmadillo.h>
#include "hf_struct_hmbdims.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

hmbdims InitDims (
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const bool hasIntercepts = TRUE
) {
  hmbdims dims;

  dims.N_S = Y_S.n_rows;
  dims.N_Sa = X_Sa.n_rows;
  dims.N_U = Z_U.n_rows;

  dims.P_X = X_S.n_cols - (hasIntercepts ? 1 : 0);
  dims.P_Z = Z_Sa.n_cols - (hasIntercepts ? 1 : 0);

  dims.df_S_X = dims.N_S - dims.P_X - (hasIntercepts ? 1 : 0);
  dims.df_Sa_Z = dims.N_Sa - dims.P_Z - (hasIntercepts ? 1 : 0);

  // Check row dimensions
  if (dims.N_S != (int)X_S.n_rows)
    Rcpp::stop("Row-dimensions of Y_S and X_S differ.");

  if (dims.N_Sa != (int)Z_Sa.n_rows)
    Rcpp::stop("Row-dimensions of X_Sa and Z_Sa differ.");

  // Check col dimensions
  if (dims.P_X != (int)X_Sa.n_cols - 1)
    Rcpp::stop("Col-dimensions of X_S and X_Sa differ.");

  if (dims.P_Z != (int)Z_U.n_cols - 1)
    Rcpp::stop("Col-dimensions of Z_Sa and Z_U differ.");

  return dims;
}


hmbdims InitDims_Z_S (
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_S,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const bool hasIntercepts = TRUE
) {
  hmbdims dims = InitDims(
    Y_S,
    X_S,
    X_Sa,
    Z_Sa,
    Z_U,
    hasIntercepts
  );

  // Check Z_S row dimensions
  if (dims.N_S != (int)Z_S.n_rows)
    Rcpp::stop("Row-dimensions of Z_S differs from X_S and Y_S.");

  // Check
  if (dims.P_Z != (int)Z_S.n_cols - 1)
    Rcpp::stop("Col-dimensions of Z_S differs from Z_Sa and Z_U.");

  // Set df_S_Z
  dims.df_S_Z = dims.N_S - dims.P_Z - (hasIntercepts ? 1 : 0);

  return dims;
}
