#include "bsa_binary_mcmc_result.h"

using namespace Rcpp;

// TODO: change the type of 'block' to an enum
void BSABinaryMCMCResult::change(const int block) {
  acc[block] += 1;
}

void BSABinaryMCMCResult::add(ThetaBinary theta_cur) {
  alphas.col(i) = theta_cur.alpha;
  beta_zs.col(i) = theta_cur.beta_z;
  gamma_zs.col(i) = theta_cur.gamma_z;
  tau_sqs.col(i) = theta_cur.tau_sq;
  gamma_xs[i] = theta_cur.gamma_x;
  beta_us[i] = theta_cur.beta_u;
  // alphas.push_back(wrap(theta_cur.alpha));
  // beta_zs.push_back(wrap(theta_cur.beta_z));
  // gamma_zs.push_back(wrap(theta_cur.gamma_z));
  // tau_sqs.push_back(wrap(theta_cur.tau_sq));
  // gamma_xs.push_back(theta_cur.gamma_x);
  // beta_us.push_back(theta_cur.beta_u);
  
  i += 1;
}

// TODO: avoid transposing
List BSABinaryMCMCResult::result() const {
  NumericVector acc(this->acc.begin(), this->acc.end());
  const std::string acc_names[] = {"alpha", "beta.z", "tau.sq",
                            "beta.u.gamma.x", "gamma.z"};
  acc.names() = std::vector<std::string>(acc_names, acc_names+5);
  return List::create(_["alpha"] = trans(alphas),
                      _["beta.z"] = trans(beta_zs),
                      _["gamma.z"] = trans(gamma_zs),
                      _["tau.sq"] = trans(tau_sqs),
                      _["gamma.x"] = gamma_xs,
                      _["beta.u"] = beta_us,
                      _["acc"] = acc);
}
