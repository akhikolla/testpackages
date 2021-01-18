#include <limits>

#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>

#include "bsalite.h"
#include "bsa_continuous_model.h"
#include "bsa_mcmc_result.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

RcppExport SEXP fitbsa(SEXP y_in,
                       SEXP x_in,
                       SEXP w_in,
                       SEXP a_in,
                       SEXP b_in,
                       SEXP k2_in,
                       SEXP el2_in,
                       SEXP nrep_in,
                       SEXP sampler_jump_in,
                       SEXP v2_in,
                       SEXP m_in,
                       SEXP mu_in,
                       SEXP rho_alpha_in,
                       SEXP sd_alpha_0_in,
                       SEXP sd_alpha_x_in) {
  RNGScope scope;               // keeping R's RNG in sync
  
  NumericVector y_r(y_in);
  arma::vec y(y_r.begin(), y_r.size(), false);

  NumericVector x_r(x_in);
  arma::vec x(x_r.begin(), x_r.size(), false);

  NumericMatrix w_r(w_in);
  arma::mat w(w_r.begin(), w_r.nrow(), w_r.ncol(), false);

  NumericVector a_r(a_in);
  arma::vec a(a_r.begin(), a_r.size(), false);

  NumericVector b_r(b_in);
  arma::vec b(b_r.begin(), b_r.size(), false);

  const double k2 = as<double>(k2_in);
  const double el2 = as<double>(el2_in);

  const int n_rep = as<int>(nrep_in);

  List sampler_jump_r(sampler_jump_in);
  std::vector<arma::vec> sampler_jump;
  sampler_jump.reserve(6);

  for (List::iterator it = sampler_jump_r.begin();
       it != sampler_jump_r.end(); ++it) {
    NumericVector block_jump_r = *it;
    arma::vec block_jump(block_jump_r.begin(), block_jump_r.size(), false);
    sampler_jump.push_back(block_jump);
    //std::cout << block_jump << std::endl;
  }

  NumericMatrix v2_r(v2_in);
  arma::mat v2(v2_r.begin(), v2_r.nrow(), v2_r.ncol(), false);

  NumericMatrix m_r(m_in);
  arma::mat m(m_r.begin(), m_r.nrow(), m_r.ncol(), false);

  NumericVector mu_r(mu_in);
  arma::vec mu(mu_r.begin(), mu_r.size(), false);

  const double rho_alpha = as<double>(rho_alpha_in);
  const double sd_alpha_0 = as<double>(sd_alpha_0_in);
  const double sd_alpha_x = as<double>(sd_alpha_x_in);
  const BSAContinuousModel f(x, y, w,
                             a, b, k2, el2,
                             v2, m, mu,
                             rho_alpha, 
                             sd_alpha_0, sd_alpha_x);
  BSAMCMCResult result = f(n_rep, sampler_jump);

  return result.result();
}
