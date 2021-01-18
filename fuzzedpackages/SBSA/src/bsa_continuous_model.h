#ifndef BSA_CONTINUOUS_MODEL_H
#define BSA_CONTINUOUS_MODEL_H

#include "bsa_mcmc_result.h"
#include "log_functor.h"

class BSAContinuousModel {
public:
  BSAContinuousModel(const arma::vec& x, 
                     const arma::vec& y, 
                     const arma::mat& w, 
                     const arma::vec& a, 
                     const arma::vec& b,
                     const double k2,
                     const double el2,
                     const arma::mat& sigma,
                     const arma::mat& m,
                     const arma::vec& mu,
                     const double rho_alpha, 
                     const double sd_alpha_0,
                     const double sd_alpha_x)
    : el2(el2),
      n(x.n_elem),
      p(w.n_cols),
      loglik(y, x, w, m, mu, el2),
      logpri(sigma, a, b, k2, 
             rho_alpha, sd_alpha_0, sd_alpha_x)
  {}

  BSAMCMCResult operator()(const int n_rep, const std::vector<arma::vec>& sampler_jump) const;
  
private:
  const double el2;

  const int n;
  const int p;

  const LogLik loglik;
  const LogPri logpri;
};
#endif
