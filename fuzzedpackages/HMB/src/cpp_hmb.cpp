//#define ARMA_64BIT_WORD 1
//#define ARMA_USE_CXX11

#include <RcppArmadillo.h>
#include "hf_struct_hmbdims.h"
#include "hf_InitDims.h"
#include "hf_MultInv.h"
#include "hf_MuVar.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
Rcpp::List cpp_hmb(
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U) {
  const hmbdims dims = InitDims(Y_S, X_S, X_Sa, Z_Sa, Z_U);

  // Inverses of predictors
  const mat Beta_inv_mat = MultInv(X_S.t() * X_S);
  const mat Alpha_inv_mat = MultInv(Z_Sa.t() * Z_Sa);

  // Storing calculation done often
  const mat G_mat_X_S = Beta_inv_mat * X_S.t();
  const mat G_mat_Z_Sa = Alpha_inv_mat * Z_Sa.t();

  // Estimated parameters
  const mat Beta = G_mat_X_S * Y_S;
  const mat Alpha = G_mat_Z_Sa * X_Sa * Beta;

  // Calculating residuals
  const vec res_S = Y_S - X_S * Beta;
  const vec res_Sa = Z_Sa * Alpha - X_Sa * Beta;
  
    // Covariances according to Eqs. (12, 34)
  const double omega2 = as_scalar(res_S.t() * res_S) / dims.df_S_X;
  const mat BetaCov = Beta_inv_mat * omega2;
 
  //const mat H_mat = Z_Sa * G_mat_Z_Sa;
  //const double Q = trace(X_Sa * BetaCov * X_Sa.t() * (eye(dims.N_Sa, dims.N_Sa) - H_mat));
  const double sigma2 = (as_scalar(res_Sa.t() * res_Sa) )/dims.df_Sa_Z;
  const mat AlphaCov = Alpha_inv_mat * sigma2 + G_mat_Z_Sa * X_Sa * Beta_inv_mat * X_Sa.t() * G_mat_Z_Sa.t() * omega2;
  
  // Getting mu estimation and mu variance estimation
  //long double muVar[2]; MuVar(muVar, Z_U, Alpha, AlphaCov);
  arma::vec one = ones(Z_U.n_rows, 1);
  arma::vec jota = one/Z_U.n_rows;
  const double muVar = as_scalar((jota.t()*Z_U*AlphaCov*Z_U.t()*jota));
  const double mu = as_scalar((jota.t()*Z_U*Alpha));

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
