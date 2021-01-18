#include "dblpm.h"

dblpm::dblpm (unsigned int T_, unsigned int N_, unsigned int M_, unsigned int L_, unsigned int D_, // GLOBAL PARS
              arma::mat edgelist_, // DATA
              arma::mat x_, arma::cube w_, arma::vec gamma_, arma::vec beta_, // LIKELIHOOD PARS
              double tauw_, double tauw0_, double taugamma_, double taugamma0_, double taubeta_, double taubeta0_, // OTHER PARS
              double taux_, double delta_, double aw_, double bw_, double agamma_, double bgamma_, double abeta_, double bbeta_, // HYPERPARS
              unsigned int niter_, unsigned int burnin_, unsigned int thin_, // MCMC SETTINGS
              double x_var_, double w_var_, double gamma_var_, double beta_var_) // PROPOSAL PARS
{
  if (debug) Rcpp::Rcout << "dblpm::dblpm has been called" << std::endl;
  T = T_;
  N = N_;
  M = M_;
  L = L_;
  D = D_;

  edgelist = edgelist_;

  x = x_;
  w = w_;
  gamma = gamma_;
  beta = beta_;

  tauw = tauw_;
  tauw0 = tauw0_;
  taugamma = taugamma_;
  taugamma0 = taugamma0_;
  taubeta = taubeta_;
  taubeta0 = taubeta0_;
  
  taux = taux_;
  delta = delta_;
  aw = aw_;
  bw = bw_;
  agamma = agamma_;
  bgamma = bgamma_;
  abeta = abeta_;
  bbeta = bbeta_;

  niter = niter_;
  burnin = burnin_;
  thin = thin_;

  x_var = x_var_;
  w_var = w_var_;
  gamma_var = gamma_var_;
  beta_var = beta_var_;

  EvaluateSumOfSquares();
  SetNoMissingData();
  FillActivity();
  FillY();
  Posterior();
  
  total_niter = burnin + niter*thin;
  x_store.set_size(niter);
  for (unsigned int iter=0; iter<niter; ++iter) x_store.at(iter).zeros(N,D);
  w_store.set_size(niter);
  for (unsigned int iter=0; iter<niter; ++iter) w_store.at(iter).zeros(M,D,T);
  gamma_store.zeros(niter,T);
  beta_store.zeros(niter,T);
  
  tauw_store.set_size(niter);
  tauw_store.fill(0.0);
  tauw0_store.set_size(niter);
  tauw0_store.fill(0.0);
  taugamma_store.set_size(niter);
  taugamma_store.fill(0.0);
  taugamma0_store.set_size(niter);
  taugamma0_store.fill(0.0);
  taubeta_store.set_size(niter);
  taubeta_store.fill(0.0);
  taubeta0_store.set_size(niter);
  taubeta0_store.fill(0.0);
  
  posterior_store.zeros(niter);
  
  if (debug) Rcpp::Rcout << "dblpm::dblpm has terminated" << std::endl;
  
}
