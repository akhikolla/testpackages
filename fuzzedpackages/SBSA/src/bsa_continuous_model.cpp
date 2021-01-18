#include "bsa_continuous_model.h"
#include "mcmc_sampler.h"

using namespace arma;

BSAMCMCResult
BSAContinuousModel::operator()(const int n_rep,
                               const std::vector<vec>& sampler_jump) const {
  Theta theta_cur(ones<vec>(2) * .1,  // alpha
                  ones<vec>(p) * .1,  // beta.z
                  ones<vec>(p) * .1,  // gamma.z
                  ones<vec>(p) * .01, // tau.sq
                  0, 0, 1);           // gamma.x, beta.u, sigma.sq
  double log_lik_cur = loglik(theta_cur);
  double log_pri_cur = logpri(theta_cur);

  const McmcReparametrizingSampler& alpha = ReparametrizeAlpha(loglik, logpri, sampler_jump.at(0));
  const McmcReparametrizingSampler& beta_z = ReparametrizeBetaZ(loglik, logpri, p, sampler_jump.at(1));
  const McmcReparametrizingSampler& sigma_sq = ReparametrizeSigmaSq(loglik, logpri, as_scalar(sampler_jump.at(2)));
  const McmcReparametrizingSampler& tau_sq = ReparametrizeTauSq(loglik, logpri, p, sampler_jump.at(3));
  const McmcReparametrizingSampler& beta_u_gamma_x = ReparametrizeBetaUGammaX(loglik, logpri, as_scalar(sampler_jump.at(4)), el2);
  const McmcReparametrizingSampler& gamma_z = ReparametrizeGammaZ(loglik, logpri, p, sampler_jump.at(5));
  const McmcReparametrizingSampler* samplers[] = {
    &alpha, &beta_z, &sigma_sq, &tau_sq, &beta_u_gamma_x, &gamma_z
  };

  BSAMCMCResult result(n_rep, p);

  for (int i = 0; i < n_rep; ++i) {
    for (int j = 0; j < 6; ++j) {
      bool changed = samplers[j]->sample(theta_cur, log_lik_cur, log_pri_cur);
      if (changed) {
        result.change(j);
      }
    }
    result.add(theta_cur);
  }

  return result;
}
