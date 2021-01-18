//#define ARMA_64BIT_WORD 1
//#define ARMA_USE_CXX11

#include <RcppArmadillo.h>
#include "cpp_gtsmb_inner.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
Rcpp::List cpp_gtsmb(
  const arma::vec& y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const arma::mat& Omega_S,
  Rcpp::NumericVector& Phi_Sa
) {
  // Put the Sa-covariances in a cube (3d-vector-ish)
  const IntegerVector Phi_Sa_dims = Phi_Sa.attr("dim");
  const cube Phi_Sa_cube(
    Phi_Sa.begin(),
    Phi_Sa_dims[0], Phi_Sa_dims[1], Phi_Sa_dims[2],
    false
  );

  // Call the inner function
  const List ret = cpp_gtsmb_inner(
    y_S,
    X_S,
    X_Sa,
    Z_Sa,
    Z_U,
    Omega_S,
    Phi_Sa_cube
  );

  return ret;
}
