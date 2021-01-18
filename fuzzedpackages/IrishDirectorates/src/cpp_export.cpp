#include "dblpm.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_dblpm_posterior(unsigned int T, unsigned int N, unsigned int M, unsigned int L, unsigned int D, // GLOBAL PARS
                               arma::mat edgelist, // DATA
                               arma::mat x, arma::cube w, arma::vec gamma, arma::vec beta, // LIKELIHOOD PARS
                               double tauw, double tauw0, double taugamma, double taugamma0, double taubeta, double taubeta0, // OTHER PARS
                               double taux, double delta, double aw, double bw, double agamma, double bgamma, double abeta, double bbeta) // HYPERPARS
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  dblpm network(T, N, M, L, D, // GLOBAL PARS
                edgelist, // DATA
                x, w, gamma, beta, // LIKELIHOOD PARS
                tauw, tauw0, taugamma, taugamma0, taubeta, taubeta0, // OTHER PARS
                taux, delta, aw, bw, agamma, bgamma, abeta, bbeta, // HYPERPARS);
                0, 0, 1, // MCMC SETTINGS
                1, 1, 1, 1); // PROPOSAL PARS
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("likelihood_value") = network.likelihood_value,
                             Rcpp::Named("posterior_value") = network.posterior_value));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_dblpm_mcmc(unsigned int T, unsigned int N, unsigned int M, unsigned int L, unsigned int D, // GLOBAL PARS
                          arma::mat edgelist, // DATA
                          arma::mat x, arma::cube w, arma::vec gamma, arma::vec beta, // LIKELIHOOD PARS
                          double tauw, double tauw0, double taugamma, double taugamma0, double taubeta, double taubeta0, // OTHER PARS
                          double taux, double delta, double aw, double bw, double agamma, double bgamma, double abeta, double bbeta, // HYPERPARS
                          unsigned int niter, unsigned int burnin, unsigned int thin, // MCMC SETTINGS
                          double x_var, double w_var, double gamma_var, double beta_var,
                          bool verbose) // PROPOSAL PARS
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  dblpm network(T, N, M, L, D, // GLOBAL PARS
                edgelist, // DATA
                x, w, gamma, beta, // LIKELIHOOD PARS
                tauw, tauw0, taugamma, taugamma0, taubeta, taubeta0, // OTHER PARS
                taux, delta, aw, bw, agamma, bgamma, abeta, bbeta, // HYPERPARS);
                niter, burnin, thin, // MCMC SETTINGS
                x_var, w_var, gamma_var, beta_var); // PROPOSAL PARS
  network.MCMC(verbose);
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("samples") = Rcpp::List::create(
                               Rcpp::Named("x_store") = network.x_store,
                               Rcpp::Named("w_store") = network.w_store,
                               Rcpp::Named("gamma_store") = network.gamma_store,
                               Rcpp::Named("beta_store") = network.beta_store,
                               Rcpp::Named("tauw_store") = network.tauw_store,
                               Rcpp::Named("tauw0_store") = network.tauw0_store,
                               Rcpp::Named("taugamma_store") = network.taugamma_store,
                               Rcpp::Named("taugamma0_store") = network.taugamma0_store,
                               Rcpp::Named("taubeta_store") = network.taubeta_store,
                               Rcpp::Named("taubeta0_store") = network.taubeta0_store),
                               // Rcpp::Named("posterior_values") = network.posterior_store,
                               Rcpp::Named("tail") = Rcpp::List::create(
                                 Rcpp::Named("dim") = Rcpp::List::create(
                                   Rcpp::Named("tframes") = network.T,
                                   Rcpp::Named("N") = network.N,
                                   Rcpp::Named("M") = network.M,
                                   Rcpp::Named("D") = network.D,
                                   Rcpp::Named("L") = network.L),
                                   Rcpp::Named("edgelist") = network.edgelist+1,
                                   Rcpp::Named("like_pars") = Rcpp::List::create(
                                     Rcpp::Named("x") = network.x,
                                     Rcpp::Named("w") = network.w,
                                     Rcpp::Named("gamma") = network.gamma,
                                     Rcpp::Named("beta") = network.beta),
                                     Rcpp::Named("tau_pars") = Rcpp::List::create(
                                       Rcpp::Named("tauw") = network.tauw,
                                       Rcpp::Named("tauw0") = network.tauw0,
                                       Rcpp::Named("taugamma") = network.taugamma,
                                       Rcpp::Named("taugamma0") = network.taugamma0,
                                       Rcpp::Named("taubeta") = network.taubeta,
                                       Rcpp::Named("taubeta0") = network.taubeta0),
                                       Rcpp::Named("hyper_pars") = Rcpp::List::create(
                                         Rcpp::Named("taux") = network.taux,
                                         Rcpp::Named("delta") = network.delta,
                                         Rcpp::Named("aw") = network.aw,
                                         Rcpp::Named("bw") = network.bw,
                                         Rcpp::Named("agamma") = network.agamma,
                                         Rcpp::Named("bgamma") = network.bgamma,
                                         Rcpp::Named("abeta") = network.abeta,
                                         Rcpp::Named("bbeta") = network.bbeta))));
}
