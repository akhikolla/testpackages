#ifndef BSA_LIK_H
#define BSA_LIK_H

#include "bsalite.h"

class LogLik {
public:
  LogLik(const arma::vec& y_, 
         const arma::vec& x_, 
         const arma::mat& w_, 
         const arma::mat& m, 
         const arma::mat& mu, 
         const double el2_)
    : y(y_), x(x_), w(w_), 
      m(m), m_inv(inv(m)), 
      mu(mu), 
      n(x_.n_elem), p(w_.n_cols),
      el2(el2_)
  {}

  double operator ()(const Theta& theta) const;
private:
  const arma::vec y;
  const arma::vec x;
  const arma::mat w;
  const arma::mat m;
  const arma::mat m_inv;
  const arma::mat mu;
  const int n; // n_elem is really an uint, but need an int for arithmetic ops
  const unsigned int p;
  const double el2;
};

class LogPri {
public:
  LogPri(const arma::mat& sigma, 
         const arma::vec& a_, const arma::vec& b_, const double k2_,
         const double rho_alpha, 
         const double sd_alpha_0,
         const double sd_alpha_x)
    : sigma(sigma), p(sigma.n_cols-1), a(a_), b(b_), k2(k2_),
      rho_alpha(rho_alpha), 
      sd_alpha_0(sd_alpha_0), sd_alpha_x(sd_alpha_x)
  {}

  double operator ()(const Theta&) const;
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
