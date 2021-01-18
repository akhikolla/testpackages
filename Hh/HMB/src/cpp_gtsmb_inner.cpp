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

Rcpp::List cpp_gtsmb_inner(
  const arma::vec& y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const arma::mat& Omega_S,
  const arma::cube& Phi_Sa
) {
  // Assuming forms:
  // y = X Beta + epsilon, epsilon ~ N(0, omega_const^2 * Omega)
  // y = Z Beta + upsilon, upsilon ~ N(0, sigma_const^2 * Sigma)
  // x_k = Z gamma_k + nu_k, nu_k ~ N(0, phi_const^2 * Phi_k)

  const hmbdims dims = InitDims(y_S, X_S, X_Sa, Z_Sa, Z_U);

  // Calculate Beta-related stuff
  mat BetaCov;
  vec Beta;
  double omega_const;
  {
    const mat Omega_S_inv = MultInv(Omega_S);
    BetaCov = MultInv(X_S.t() * Omega_S_inv * X_S);
    Beta = BetaCov * X_S.t() * Omega_S_inv * y_S;

    // Calculate omega_const
    const vec resid_S_X = y_S - X_S * Beta;
    omega_const = as_scalar(resid_S_X.t() * Omega_S_inv * resid_S_X) / dims.df_S_X;

    // Fix BetaCov
    BetaCov = BetaCov * omega_const;
  }

  // Declare some stuff for use in loop. Also Declare GammaCov_ish.
  // GammaCov_ish corresponds to the joint terms of the first (1) and third (3)
  // term of Eq. (12) in Holm et. al. (2017).
  cube Phi_Sa_inv(dims.N_Sa, dims.N_Sa, dims.P_X);
  mat phi_const_mat = zeros(dims.P_X, dims.P_X);

  mat GammaCov_ish(dims.P_Z + 1, dims.P_Z + 1, fill::zeros);

  // Declaring the Gamma matrix, also setting the first column to [1 0 0 ...]'
  mat Gamma(dims.P_Z + 1, dims.P_X + 1, fill::zeros);
  Gamma.at(0, 0) = 1;

  // Loop through all the regressors (except intercept)
  for (int p = 1; p <= dims.P_X; ++p) {
    int rowstart = p * (p - 1) / 2;
    // Store inverses for easy access (used here and in inner loop).
    // Also store some calculations for easy access.
    Phi_Sa_inv.slice(p - 1) = MultInv(Phi_Sa.slice(rowstart + p - 1));

    const mat G_Z_left = MultInv(Z_Sa.t() * Phi_Sa_inv.slice(p - 1) * Z_Sa)
      * Z_Sa.t();

    // Calculate the current Gamma column.
    Gamma.col(p) = G_Z_left * Phi_Sa_inv.slice(p - 1) * X_Sa.col(p);

    // Calculate the GammaCov_ish for the (p, p)-case
    // Comparing with Holm et. al. (2017) Eq. (12), the lower row of the
    // expression corresponds to the GammaCov-elements
    // (omega-covs in Holm et. al.)
    const vec resid_X_Z_p = X_Sa.col(p) - Z_Sa * Gamma.col(p);
    const double phi_const_pp = as_scalar(
      resid_X_Z_p.t() * Phi_Sa_inv.slice(p - 1) * resid_X_Z_p
    ) / dims.df_Sa_Z;
    phi_const_mat.at(p - 1, p - 1) = phi_const_pp;

    GammaCov_ish += phi_const_pp
        * (Beta.at(p) * Beta.at(p) - BetaCov.at(p, p))
        * G_Z_left * Phi_Sa_inv.slice(p - 1) * G_Z_left.t();

    // Loop through all the regressors below the current one. Since this loop
    // and other calculations are symmetric, we don't have to loop through all
    // regressors in this step, just multiply the results by two.
    for (int k = 1; k < p; ++k) {
      const mat interPhi = Phi_Sa.slice(rowstart + k - 1);

      // Calculating phi_const for p/k.
      const vec resid_X_Z_k = X_Sa.col(k) - Z_Sa * Gamma.col(k);
      const double phi_const_pk = as_scalar(
        resid_X_Z_p.t() * MultInv(interPhi) * resid_X_Z_k
      ) / dims.df_Sa_Z;
      phi_const_mat.at(p - 1, k - 1) = phi_const_pk;

      const mat G_Z_right = Z_Sa
        * MultInv(Z_Sa.t() * Phi_Sa_inv.slice(k - 1) * Z_Sa);

      // Doubled due to symmetry of involved elements.
      GammaCov_ish += phi_const_pk
        * 2 * (Beta.at(p) * Beta.at(k))
        * G_Z_left
        * Phi_Sa_inv.slice(p - 1)
        * interPhi
        * Phi_Sa_inv.slice(k - 1)
        * G_Z_right;
        /*GammaCov_ish += phi_const_pk
        * 2 * (Beta.at(p) * Beta.at(k))
        * G_Z_left
        * Str_SqMult(Str_SqMult(Phi_Sa_inv.slice(p - 1), interPhi), Phi_Sa_inv.slice(k - 1))
        * G_Z_right;*/
    }
  }

  // Calculating mu-estimator and variance-estimator of mu-estimator
  // Variance according to Holm et. al (2017): 2 + (1 - 3)
  //long double muVar[2];
  //MuVar(muVar, Z_U, Gamma * Beta, (Gamma * BetaCov * Gamma.t() + GammaCov_ish));
  
  arma::vec one = ones(Z_U.n_rows, 1);
  arma::vec jota = one/Z_U.n_rows;
  const double muVar = as_scalar((jota.t()*Z_U*(Gamma * BetaCov * Gamma.t() + GammaCov_ish)*Z_U.t()*jota));
  const double mu = as_scalar((jota.t()*Z_U*Gamma*Beta));

  List ret;
  ret["Beta"] = Beta;
  ret["BetaCov"] = BetaCov;
  ret["omega2"] = omega_const;
  ret["Gamma"] = Gamma;
  ret["mu"] = mu;
  ret["muVar"] = muVar;
  ret["GammaCov_ish"] = GammaCov_ish;
  ret["phi2s"] = phi_const_mat;
  return ret;
}
