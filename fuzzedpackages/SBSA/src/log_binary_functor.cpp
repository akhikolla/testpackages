#include <limits>

#include "log_binary_functor.h"

#include <R.h>
#include <Rmath.h>

using namespace arma;

double LogBinaryLik::operator ()(const ThetaBinary& theta) const {
  const double alpha_x_star = theta.alpha(1) + 
    theta.beta_u * theta.gamma_x;

  const vec beta_z_star = theta.beta_z + 
    theta.beta_u * theta.gamma_z;

  const mat D = diagmat(theta.tau_sq);
  if (min(eig_sym(m-D)) <= 1E-4) {
    // (M-D) is ill-conditioned, nothing to do here...
    return std::numeric_limits<double>::quiet_NaN();
  }
  const mat MD_inv = inv(m - D);
  const mat D_inv = inv(D);
  const mat MDD_inv = inv(MD_inv + D_inv);
  
  const double omega_0 = theta.alpha(0);

  const double omega_x = as_scalar(alpha_x_star +
                                   trans(beta_z_star) * MDD_inv * MD_inv * mu);
  const rowvec omega_w = trans(beta_z_star) * MDD_inv * D_inv;

  const double eta_sq = (theta.beta_u * theta.beta_u * el2) +
    as_scalar(trans(beta_z_star) * MDD_inv * beta_z_star);
  const double eta = sqrt(eta_sq);
  
  double llik = 0;
  for (unsigned int i = 0; i < n; ++i) {
    const double delta = as_scalar(omega_0 + omega_x * x(i) +
				   omega_w * trans(w.row(i)));
    double q = q_val(delta, eta);
    
    llik += y(i) ? log(q) : log(1-q);
  }
  
  return llik;
}

double LogBinaryLik::q_val(const double delta, const double eta) const {
  double q = 0;
  for (int i = 0; i < q_steps; ++i) {
    q += 1 / (1 + exp(-(delta + eta * q_norm(i))));
  }
  return q / q_steps;
}

static double square(const double x) {
  return x * x;
}

static vec dbeta(const vec& x, const vec& a, const vec& b) {
  vec y(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; ++i) {
    y(i) = ::Rf_dbeta(x(i), a(i), b(i), 0);
  }
  return y;
}

// Unchanged from LogPri
double LogBinaryPri::operator ()(const ThetaBinary& theta) const {
  vec beta_fk = zeros<vec>(p+1);
  beta_fk.rows(0,p-1) = theta.beta_z;
  beta_fk.row(p) = theta.beta_u;

  // beta contribution
  double beta_df = 10;
  double beta_sc = log(5.0);
  double log_pri = -(beta_df + p + 1) / 2 * 
    log(1+(1/beta_df) * accu(beta_fk%beta_fk) / beta_sc);

  // alpha contribution
  const double alpha_0 = theta.alpha(0);
  const double alpha_x = theta.alpha(1);
  log_pri += -1/(1-square(rho_alpha)) *
    (square(alpha_0/sd_alpha_0) + 
     square(alpha_x/sd_alpha_x) -
     2 * rho_alpha * alpha_0 * alpha_x / (sd_alpha_0*sd_alpha_x));

  // tau_sq contribution
  log_pri += accu(log(dbeta(theta.tau_sq, a, b)));

  vec gamma_fk = zeros<vec>(p+1);
  gamma_fk(0) = theta.gamma_x;
  gamma_fk.rows(1, p) = theta.gamma_z;

  // gamma contribution
  vec zero_tau_sq_diag = zeros<vec>(p+1);
  zero_tau_sq_diag.rows(1,p) = theta.tau_sq;
  mat tmp = (sigma - diagmat(zero_tau_sq_diag)) / k2;
  // TODO: check that tmp is always symmetric, otherwise use eig_gen
  if (min(eig_sym(tmp)) <= 1E-4) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  log_pri += accu(log(mat(chol(tmp)).diag())) - 
    0.5 * as_scalar(trans(gamma_fk) * tmp * gamma_fk);

  return log_pri;
}
