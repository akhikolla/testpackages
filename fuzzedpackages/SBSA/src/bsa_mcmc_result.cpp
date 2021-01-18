#include "bsa_mcmc_result.h"

using namespace Rcpp;

// TODO: change the type of 'block' to an enum
void BSAMCMCResult::change(const int block) {
  acc[block] += 1;
}

void BSAMCMCResult::add(Theta theta_cur) {
  alphas.col(i) = theta_cur.alpha;
  beta_zs.col(i) = theta_cur.beta_z;
  gamma_zs.col(i) = theta_cur.gamma_z;
  tau_sqs.col(i) = theta_cur.tau_sq;
  gamma_xs[i] = theta_cur.gamma_x;
  beta_us[i] = theta_cur.beta_u;
  sigma_sqs[i] = theta_cur.sigma_sq;
  // alphas.push_back(wrap(theta_cur.alpha));
  // beta_zs.push_back(wrap(theta_cur.beta_z));
  // gamma_zs.push_back(wrap(theta_cur.gamma_z));
  // tau_sqs.push_back(wrap(theta_cur.tau_sq));
  // gamma_xs.push_back(theta_cur.gamma_x);
  // beta_us.push_back(theta_cur.beta_u);
  // sigma_sqs.push_back(theta_cur.sigma_sq);
  
  i += 1;
}

// TODO: avoid transposing
List BSAMCMCResult::result() const {
  NumericVector acc(this->acc.begin(), this->acc.end());
  const std::string acc_names[] = {"alpha", "beta.z", "sigma.sq", "tau.sq",
                            "beta.u.gamma.x", "gamma.z"};
  acc.names() = std::vector<std::string>(acc_names, acc_names+6);
  return List::create(_["alpha"] = trans(alphas),
                      _["beta.z"] = trans(beta_zs),
                      _["gamma.z"] = trans(gamma_zs),
                      _["tau.sq"] = trans(tau_sqs),
                      _["gamma.x"] = gamma_xs,
                      _["beta.u"] = beta_us,
                      _["sigma.sq"] = sigma_sqs,
                      _["acc"] = acc);
}
