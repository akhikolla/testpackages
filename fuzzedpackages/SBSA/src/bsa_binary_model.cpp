#include "bsa_binary_model.h"
#include "mcmc_binary_sampler.h"

using namespace arma;

BSABinaryMCMCResult
BSABinaryModel::operator()(const int n_rep,
                           const std::vector<vec>& sampler_jump) const {
  ThetaBinary theta_cur(ones<vec>(2) * .1, // alpha
                        ones<vec>(p) * .1, // beta.z
                        ones<vec>(p) * .1, // gamma.z
                        ones<vec>(p) * .01, // tau.sq
                        0, 0);              // gamma.x, beta.u
  double log_lik_cur = loglik(theta_cur);
  double log_pri_cur = logpri(theta_cur);

  const McmcReparametrizingBinarySampler& alpha = ReparametrizeBinaryAlpha(loglik, logpri, sampler_jump.at(0));
  const McmcReparametrizingBinarySampler& beta_z = ReparametrizeBinaryBetaZ(loglik, logpri, p, sampler_jump.at(1));
  const McmcReparametrizingBinarySampler& tau_sq = ReparametrizeBinaryTauSq(loglik, logpri, p, sampler_jump.at(2));
  const McmcReparametrizingBinarySampler& beta_u_gamma_x = ReparametrizeBinaryBetaUGammaX(loglik, logpri, as_scalar(sampler_jump.at(3)), el2);
  const McmcReparametrizingBinarySampler& gamma_z = ReparametrizeBinaryGammaZ(loglik, logpri, p, sampler_jump.at(4));
  const McmcReparametrizingBinarySampler* samplers[] = {
    &alpha, &beta_z, &tau_sq, &beta_u_gamma_x, &gamma_z
  };

  BSABinaryMCMCResult result(n_rep, p);

  for (int i = 0; i < n_rep; ++i) {
    for (int j = 0; j < 5; ++j) {
      bool changed = samplers[j]->sample(theta_cur, log_lik_cur, log_pri_cur);
      if (changed) {
        result.change(j);
      }
    }
    result.add(theta_cur);
  }

  return result;
}
