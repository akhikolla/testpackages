#include "expsbm.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_expSBM_ELBO(unsigned int N, arma::mat edgelist, arma::mat Z, arma::vec lambda, arma::mat mu, arma::mat nu, bool directed, bool trunc, bool verbose)
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  expsbm network(N, edgelist, Z, lambda, mu, nu, directed, trunc, 0.001, 1, verbose);
  if (network.verbose) network.Print();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("elbo_value") = network.elbo_value));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_expSBM_EM(unsigned int N, arma::mat edgelist, arma::mat Z, arma::vec lambda, arma::mat mu, arma::mat nu, bool directed, bool trunc, double tol, unsigned int n_iter_max, bool verbose)
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  expsbm network(N, edgelist, Z, lambda, mu, nu, directed, trunc, tol, n_iter_max, verbose);
  network.Optimisation();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("elbo_values") = network.elbo_values_store,
                             Rcpp::Named("Z_star") = network.Z,
                             Rcpp::Named("lambda_star") = network.lambda,
                             Rcpp::Named("mu_star") = network.mu,
                             Rcpp::Named("nu_star") = network.nu));
}
