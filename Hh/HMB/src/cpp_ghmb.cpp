//#define ARMA_64BIT_WORD 1
//#define ARMA_USE_CXX11

#include <RcppArmadillo.h>
#include "hf_struct_hmbdims.h"
#include "hf_InitDims.h"
#include "hf_MuVar.h"
#include "hf_MultInv.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
Rcpp::List cpp_ghmb(
  const arma::vec& y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const arma::mat& Omega_S,
  const arma::mat& Sigma_Sa) {

  const hmbdims dims = InitDims(y_S, X_S, X_Sa, Z_Sa, Z_U);

  // Store for easy access
  const mat Omega_S_inv = MultInv(Omega_S);
  const mat Sigma_Sa_inv = MultInv(Sigma_Sa);
  const mat Alpha_inv_mat = MultInv(Z_Sa.t() * Sigma_Sa_inv * Z_Sa);
  const mat Beta_inv_mat = MultInv(X_S.t() * Omega_S_inv * X_S);
  const mat G_mat_X = Beta_inv_mat * X_S.t() * Omega_S_inv;
  const mat G_mat_Z = Alpha_inv_mat * Z_Sa.t() * Sigma_Sa_inv;

  // Prepare very useful matrices


  // Calculate Beta and Alpha
  const vec Beta = G_mat_X * y_S;
  const vec Alpha = G_mat_Z * X_Sa * Beta;

  // Calculate variance-constants
  double omega2, sigma2, Q;

  // Calculate omega_const
  {
    const vec resid_S = y_S - X_S * Beta;
    omega2 = as_scalar(resid_S.t() * Omega_S_inv * resid_S) / dims.df_S_X;
  }

  // Calculate BetaCov
  const mat BetaCov = omega2 * Beta_inv_mat;
  
  // Calculate sigma_const
  mat BD_mat(dims.N_Sa, dims.N_Sa, fill::zeros);
  
  {
    //const mat H_mat = Z_Sa * G_mat_Z;
    //Q = trace(X_Sa * BetaCov * X_Sa.t() * Sigma_Sa_inv * (eye(dims.N_Sa, dims.N_Sa) - H_mat));
    const vec resid_Sa = X_Sa * Beta - Z_Sa * Alpha;
    sigma2 = (as_scalar(resid_Sa.t() * Sigma_Sa_inv * resid_Sa))/dims.df_Sa_Z;
  }

  // Calculate  AlphaCov
  const mat AlphaCov = Alpha_inv_mat * sigma2 + G_mat_Z * X_Sa * Beta_inv_mat * X_Sa.t() * G_mat_Z.t() * omega2;

  // Create array for mu-estimator and variance-estimator of mu-estimator
  //long double muVar[2]; MuVar(muVar, Z_U, Alpha, AlphaCov);
  arma::vec one = ones(Z_U.n_rows, 1);
  arma::vec jota = one/Z_U.n_rows;
  const double muVar = as_scalar((jota.t()*Z_U*AlphaCov*Z_U.t()*jota));
  const double mu = as_scalar((jota.t()*Z_U*Alpha));

  // Return results
  List ret;
  ret["Beta"] = Beta;
  ret["Alpha"] = Alpha;
  ret["BetaCov"] = BetaCov;
  ret["AlphaCov"] = AlphaCov;
  ret["omega2"] = omega2;
  ret["sigma2"] = sigma2;
  ret["mu"] = mu;
  ret["muVar"] = muVar;
  return ret;
}
