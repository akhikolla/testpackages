#include "mcmc_sampler.h"

#include <R.h>
#include <Rmath.h>

using namespace arma;

bool
McmcReparametrizingSampler::sample(Theta& theta_cur, 
                                   double& log_lik_cur, 
                                   double& log_pri_cur) const {
  Theta theta_cnd = reparametrize(theta_cur);
  if (theta_cnd.gamma_x > 0) {
    double log_lik_cnd = loglik(theta_cnd);
    double log_pri_cnd = logpri(theta_cnd);

    if (!ISNAN(log_lik_cnd) && 
        min(theta_cnd.tau_sq) > 0 &&
        theta_cnd.sigma_sq > 0 &&
        log(::Rf_runif(0, 1)) < (log_lik_cnd + log_pri_cnd - 
                                 log_lik_cur - log_pri_cur)) {
      theta_cur = theta_cnd;
      log_lik_cur = log_lik_cnd;
      log_pri_cur = log_pri_cnd;
      return true;
    }
  }
  return false;
}


Theta ReparametrizeAlpha::reparametrize(const Theta& theta_cur) const {
  Theta theta(theta_cur);
  theta.alpha(0) += ::Rf_rnorm(0, sigma(0));
  theta.alpha(1) += ::Rf_rnorm(0, sigma(1));

  return theta;
}


Theta ReparametrizeBetaZ::reparametrize(const Theta& theta_cur) const {
  Theta theta(theta_cur);
  for (int i = 0; i < p; ++i) {
    theta.beta_z(i) += ::Rf_rnorm(0, sigma(i));
  }

  return theta;
}


Theta ReparametrizeSigmaSq::reparametrize(const Theta& theta_cur) const {
  Theta theta(theta_cur);
  theta.sigma_sq += ::Rf_rnorm(0, sigma);

  return theta;
}


Theta ReparametrizeTauSq::reparametrize(const Theta& theta_cur) const {
  Theta theta(theta_cur);
  for (int i = 0; i < p; ++i) {
    theta.tau_sq(i) += ::Rf_rnorm(0, sigma(i));
  }

  return theta;
}


Theta ReparametrizeBetaUGammaX::reparametrize(const Theta& theta_cur) const {
  const double beta_u_noise = ::Rf_rnorm(0, sigma);
  const double gamma_u_noise = ::Rf_rnorm(0, sigma);

  Theta theta(theta_cur);
  theta.alpha(1) -= theta.beta_u * gamma_u_noise +
    theta.gamma_x * beta_u_noise +
    beta_u_noise * gamma_u_noise;
  theta.beta_z -= beta_u_noise * theta.gamma_z;
  theta.sigma_sq -= el2 * (2 * beta_u_noise * theta.beta_u +
                           beta_u_noise * beta_u_noise);
  theta.beta_u += beta_u_noise;
  theta.gamma_x += gamma_u_noise;

  return theta;
}


Theta ReparametrizeGammaZ::reparametrize(const Theta& theta_cur) const {
  vec gamma_z_noise = zeros<vec>(p);
  for (int i = 0; i < p; ++i) {
    gamma_z_noise(i) = ::Rf_rnorm(0, sigma(i));
  }

  Theta theta(theta_cur);
  theta.gamma_z += gamma_z_noise;
  theta.beta_z -= theta_cur.beta_u * gamma_z_noise;
  return theta;
}
