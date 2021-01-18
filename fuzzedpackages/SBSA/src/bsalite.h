#ifndef BSA_LITE_H
#define BSA_LITE_H

#include <RcppArmadillo.h>

struct Theta {
  arma::vec alpha;
  arma::vec beta_z;
  arma::vec gamma_z;
  arma::vec tau_sq;
  double gamma_x;
  double beta_u;
  double sigma_sq;

  Theta(arma::vec alpha_, arma::vec beta_z_, 
	arma::vec gamma_z_, arma::vec tau_sq_,
	double gamma_x_, double beta_u_, double sigma_sq_)
    : alpha(alpha_),
      beta_z(beta_z_),
      gamma_z(gamma_z_),
      tau_sq(tau_sq_),
      gamma_x(gamma_x_),
      beta_u(beta_u_),
      sigma_sq(sigma_sq_)
  {}
};

struct ThetaBinary {
  arma::vec alpha;
  arma::vec beta_z;
  arma::vec gamma_z;
  arma::vec tau_sq;
  double gamma_x;
  double beta_u;

  ThetaBinary(arma::vec alpha_, arma::vec beta_z_, 
	      arma::vec gamma_z_, arma::vec tau_sq_,
	      double gamma_x_, double beta_u_)
    : alpha(alpha_),
      beta_z(beta_z_),
      gamma_z(gamma_z_),
      tau_sq(tau_sq_),
      gamma_x(gamma_x_),
      beta_u(beta_u_)
  {}
};

#endif
