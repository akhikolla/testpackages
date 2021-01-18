#include <limits>

#include "log_functor.h"

#include <R.h>
#include <Rmath.h>

using namespace arma;

double LogLik::operator ()(const Theta& theta) const {
  mat D = diagmat(theta.tau_sq);
  
  // TODO: is V-D for sure symmetric? if not, use eig_gen instead of eig_sym
  if (min(eig_sym(m-D)) <= 1E-4) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double alpha_x_star = theta.alpha(1) + 
    theta.beta_u * theta.gamma_x;

  const vec beta_z_star = theta.beta_z + 
    theta.beta_u * theta.gamma_z;

  vec coef = zeros(2+p);
  coef(0) = theta.alpha(0);
  coef.row(1) = as_scalar(alpha_x_star + 
                           trans(beta_z_star) * D * m_inv * mu);
  coef.rows(2, coef.n_elem-1) = trans(trans(beta_z_star) * (eye(p, p) - D * m_inv));

  const double sigma_sq_star = theta.sigma_sq + 
    theta.beta_u * theta.beta_u * el2;
  double vr = as_scalar(sigma_sq_star + 
                        trans(beta_z_star) * 
                        (D - D*m_inv*D) * 
                        beta_z_star);

  mat one_x_w = ones(x.n_rows, p + 2);
  one_x_w.col(1) = x;
  one_x_w.cols(2, p+1) = w;
  vec foo = y - one_x_w * coef;
  double llik = as_scalar((-n/2) * log(vr) - sum(square(foo)) / (vr * 2));
  return llik;
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

double LogPri::operator ()(const Theta& theta) const {
  vec beta_fk = zeros<vec>(p+1);
  beta_fk.rows(0,p-1) = theta.beta_z;
  beta_fk.row(p) = theta.beta_u;

  // beta contribution
  double beta_df = 10;
  double beta_sc = log(5.0);
  double log_pri = -(beta_df + p + 1) / 2 * 
    log(1+(1/beta_df) * accu(square(beta_fk)) / beta_sc);

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
