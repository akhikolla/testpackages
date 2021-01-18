#ifndef BSA_BINARY_MODEL_H
#define BSA_BINARY_MODEL_H

#include "bsa_binary_mcmc_result.h"
#include "log_binary_functor.h"

class BSABinaryModel {
public:
  BSABinaryModel(const arma::vec& x, 
                 const arma::vec& y, 
                 const arma::mat& w, 
                 const arma::vec& a, 
                 const arma::vec& b,
                 const double k2,
                 const double el2,
                 const arma::mat& sigma,
                 const arma::mat& m,
                 const arma::vec& mu,
                 const double q_steps,
                 const arma::vec& q_norm,
                 const double rho_alpha, 
                 const double sd_alpha_0,
                 const double sd_alpha_x)
    : el2(el2),
      n(x.n_elem),
      p(w.n_cols),
      loglik(y, x, w, m, mu, el2, q_steps, q_norm),
      logpri(sigma, a, b, k2, 
             rho_alpha, sd_alpha_0, sd_alpha_x)
  {}
  
  BSABinaryMCMCResult operator()(const int n_rep, const std::vector<arma::vec>& sampler_jump) const;
  
private:
  const double el2;

  const int n;
  const int p;

  const LogBinaryLik loglik;
  const LogBinaryPri logpri;
};
#endif
