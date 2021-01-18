#ifndef BSA_LIK_BIN_H
#define BSA_LIK_BIN_H

#include "bsalite.h"

class LogBinaryLik {
public:
  LogBinaryLik(const arma::vec& y_,
               const arma::vec& x_,
               const arma::mat& w_,
               const arma::mat& m,
               const arma::mat& mu,
               const double el2_,
               const double q_steps,
               const arma::vec& q_norm)
    : y(y_), x(x_), w(w_), 
      m(m), mu(mu), 
      n(x_.n_elem), p(w_.n_cols),
      el2(el2_),
      q_steps(q_steps), q_norm(q_norm)
  {}

  double operator ()(const ThetaBinary& theta) const;
private:
  double q_val(const double delta, const double eta) const;
  
  const arma::vec y;
  const arma::vec x;
  const arma::mat w;
  const arma::mat m;
  const arma::mat mu;
  const unsigned int n;
  const unsigned int p;
  const double el2;

  const int q_steps;      // number of samples for integration
  const arma::vec q_norm; // quantiles of normal := qnorm((1:q_steps)/(q_steps+1))
};

class LogBinaryPri {
public:
  LogBinaryPri(const arma::mat& sigma, 
               const arma::vec& a_, const arma::vec& b_, const double k2_,
               const double rho_alpha,
               const double sd_alpha_0,
               const double sd_alpha_x)
    : sigma(sigma), p(sigma.n_cols-1), a(a_), b(b_), k2(k2_),
      rho_alpha(rho_alpha), 
      sd_alpha_0(sd_alpha_0), sd_alpha_x(sd_alpha_x)
  {}

  double operator ()(const ThetaBinary&) const;
private:
  const arma::mat sigma;
  const unsigned int p;
  const arma::vec a; 
  const arma::vec b;
  const double k2;

  // Parameters for the prior of alpha
  const double rho_alpha;
  const double sd_alpha_0;
  const double sd_alpha_x;
};

#endif
