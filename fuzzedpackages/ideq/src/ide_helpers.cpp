#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "distributions.h"

using namespace Rcpp;

double kernelLikelihood(const arma::mat & G, const arma::mat & theta,
                        const arma::mat & W) {
  const int T = theta.n_cols-1;
  arma::mat tmp = arma::zeros(1, 1);
  arma::colvec d;
  arma::mat W_inv = arma::inv_sympd(W);
  
  for (int t=1; t<=T; ++t) {
    d = theta.col(t) - G * theta.col(t-1);
    tmp += d.t() * W_inv * d;
  }
  
  return -arma::as_scalar(tmp)/2.0;
}

double kernelLikelihoodDiscount(const arma::mat & G, const arma::mat & theta, 
                        const arma::cube & C, const double lambda) {
  const int T = theta.n_cols-1;
  arma::mat tmp = arma::zeros(1, 1);
  arma::colvec d;
  
  for (int t=1; t<=T; ++t) {
    d = theta.col(t) - G * theta.col(t-1);
    tmp += d.t() * arma::solve(G * C.slice(t) * G.t(), d) / lambda;
  }
  
  return -arma::as_scalar(tmp)/2.0;
}

arma::mat makeW(const int J, const double L) {
  // NOTE: function is assumed to have period of L
  arma::colvec freqs1 = 2*PI/L * arma::regspace(1, J);
  arma::colvec freqs2 = 2*PI/L * arma::regspace(-J,J);
  arma::mat w(2*J*(J+1), 2);
  
  // Create w
  w.col(0).rows(0, (2*J+1)*J-1) = arma::repmat(freqs2, J, 1);
  w.col(0).rows((2*J+1)*J, w.n_rows-1) = freqs1;
  
  w.col(1).rows(0, (2*J+1)*J-1) = arma::repelem(freqs1, 2*J+1, 1);
  w.col(1).rows((2*J+1)*J, w.n_rows-1) = arma::zeros(J, 1);
  
  return w;
}

// Observation matrix
// The number of total basis function is J^2+1
// L is the range of the Fourier approximation
// locs are the centered/scaled spatial locations
arma::mat makeF(const arma::mat & locs, const arma::mat & w, const double L) {
  arma::mat Jmat = locs.col(0) * w.col(0).t() +
                   locs.col(1) * w.col(1).t();
  const int k = Jmat.n_cols;
  arma::mat Phi(Jmat.n_rows, 2*k + 1);
  Phi.col(0).fill(1);
  Phi.cols(1, k) = sqrt(2.0) * arma::cos(Jmat);
  Phi.cols(k + 1, 2*k) = sqrt(2.0) * arma::sin(Jmat);
  Phi /= L;
  return Phi;
}

// The function makeB returns the matrix B used as part of the process matrix
// mu and Sigma are the parameters of the IDE kernel
void makeB(arma::mat & B, const arma::mat & mu, const arma::cube & Sigma, 
           const arma::mat & locs, const arma::mat & w, const double L) {
  
  const bool SV_mu = mu.n_cols > 1;
  const bool SV_Sigma = Sigma.n_slices > 1;
  
  // Jmat1 and Jmat2
  arma::mat Jmat1, Jmat2;
  
  if (SV_mu) {
    Jmat1 = (locs.col(0) - mu.col(0)) * w.col(0).t() +
            (locs.col(1) - mu.col(1)) * w.col(1).t();
  } else {
    Jmat1 = (locs.col(0) - mu(0)) * w.col(0).t() +
            (locs.col(1) - mu(1)) * w.col(1).t();
  }
  
  if (SV_Sigma) {
    arma::colvec tmp = Sigma.tube(0, 0);
    Jmat2 = tmp * arma::square(w.col(0).t());
    tmp = Sigma.tube(1, 1);
    Jmat2 += tmp * arma::square(w.col(1).t());
    tmp = 2 * Sigma.tube(0, 1);
    Jmat2 += tmp * arma::prod(w.t(), 0);
  } else {
    arma::mat Jvec = Sigma.at(0, 0, 0) * arma::square(w.col(0)) +
                     Sigma.at(1, 1, 0) * arma::square(w.col(1)) +
                     2 * Sigma.at(0, 1, 0) * arma::prod(w, 1);
    Jmat2 = arma::kron( arma::ones(locs.n_rows, 1), Jvec.t() );
  }
  
  // Exponentiate Jmat2
  Jmat2 = arma::exp(-0.5 * Jmat2);
  
  // B
  const int k = Jmat1.n_cols;
  B.col(0).fill(1);
  B.cols(1, k) = sqrt(2.0) * Jmat2 % arma::cos(Jmat1);
  B.cols(k+1, 2*k) = sqrt(2.0) * Jmat2 % arma::sin(Jmat1);
  B /= L;
  return;
}

void mapSigma(arma::cube & s_many, const arma::cube & s_few,
              const arma::mat K) {
  if (s_many.n_rows != s_few.n_rows || s_many.n_cols != s_few.n_cols) {
    throw std::invalid_argument("s_many and s_few must have same number of rows and columns");
  }
  arma::colvec tmp;
  
  for (unsigned int r=0; r<s_few.n_rows; ++r) {
    for (unsigned int c=0; c<s_few.n_cols; ++c) {
      tmp = s_few.tube(r, c);
      s_many.tube(r, c) = K * tmp;
    }
  }
  
  return;
}

arma::mat proposeMu(arma::mat mu, arma::mat Sigma) {
  arma::colvec mu_vec = arma::vectorise(mu);
  if (!Sigma.is_square()) {
    throw std::invalid_argument("Sigma must be square");
  }
  if (mu_vec.n_elem != Sigma.n_rows) {
    throw std::invalid_argument("Number of elements in mu must equal dimension of Sigma");
  }
  
  arma::colvec tmp = rmvnorm(mu_vec, Sigma);
  arma::mat mu_proposal = arma::reshape(tmp, mu.n_rows, mu.n_cols);
  return mu_proposal;
}
