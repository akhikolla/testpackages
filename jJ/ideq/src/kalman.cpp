#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "misc_helpers.h"

using namespace Rcpp;

void kalman(arma::mat & m, arma::cube & C, arma::mat & a, arma::cube & R_inv,
            const arma::mat & Y, const arma::mat & F, const arma::mat & G,
            const double sigma2, const double lambda, const arma::mat & W) {
  // This function assumes that V is (sigma2 * I)
  const int T = Y.n_cols-1;
  const int S = Y.n_rows;
  const int P = G.n_rows;
  const bool Discount = lambda > 0;

  // Don't need to keep these quantities
  arma::mat Q(S, S), Q_inv(S, S), R_t(P, P), FR(S, P), RF_t(P, S);
  arma::colvec f(S);
  // Could also save transposed copies of F, G

  for (int t = 1; t <= T; ++t) {
    checkUserInterrupt();
    // One step ahead predictive distribution of theta
    a.col(t) = G * m.col(t-1);
    if (Discount) R_t = (1 + lambda) * G * C.slice(t-1) * G.t();
    else R_t = G * C.slice(t-1) * G.t() + W;

    // One step ahead predictive distribution of Y_t
    f = F * a.col(t);
    FR = F * R_t;
    Q = FR * F.t();
    Q.diag() += sigma2;
    
    try {
      Q_inv = arma::inv_sympd(Q);
    }
    catch (std::runtime_error & e) {
      Rcout << "Failed to invert Q in kalman filter" << std::endl;
      Rcout << "Consider making sigma2 larger" << std::endl;
      Q_inv = forceInv(Q);
    }

    // Filtering distribution of theta
    RF_t = FR.t();
    m.col(t) = a.col(t) + RF_t * Q_inv * (Y.col(t) - f);
    C.slice(t) = R_t - RF_t * Q_inv * FR;

    // Invert R for sampling
    try {
      R_inv.slice(t) = arma::inv_sympd(R_t);
    }
    catch (std::runtime_error & e) {
      Rcout << "Failed to invert R in kalman filter" << std::endl;
      Rcout << "Consider using a stronger prior for W" << std::endl;
      R_inv.slice(t) = forceInv(R_t);
    }
    R_inv.slice(t) = arma::inv_sympd(R_t);
  }
  return;
}
