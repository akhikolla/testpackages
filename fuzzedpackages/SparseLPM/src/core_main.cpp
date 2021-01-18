#include <RcppArmadillo.h>
#include "core_slpm_var.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_SLPM_ELBO(arma::mat adj, 
                         arma::mat var_alpha_u, arma::mat var_alpha_v, 
                         arma::mat var_beta_u, arma::mat var_beta_v, 
                         arma::cube var_lambda, 
                         arma::vec var_delta, arma::vec var_a, arma::vec var_b, 
                         arma::vec delta, arma::vec a_gamma, arma::vec b_gamma, 
                         bool verbose)
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  slpm_var network(adj, var_alpha_u, var_alpha_v, var_beta_u, var_beta_v, var_lambda, var_delta, var_a, var_b, delta, a_gamma, b_gamma, verbose);
  if (network.verbose) network.Print();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("elbo") = network.elbo_value));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_SLPM_Optimisation(arma::mat adj, 
                                 arma::mat var_alpha_u, arma::mat var_alpha_v, 
                                 arma::mat var_beta_u, arma::mat var_beta_v, 
                                 arma::cube var_lambda, 
                                 arma::vec var_delta, arma::vec var_a, arma::vec var_b, 
                                 arma::vec delta, arma::vec a_gamma, arma::vec b_gamma, 
                                 double tol, unsigned int n_iter_max, bool natural_gradient, 
                                 double learning_rate_factor_up, double learning_rate_factor_down, bool verbose)
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  slpm_var network(adj, var_alpha_u, var_alpha_v, var_beta_u, var_beta_v, var_lambda, var_delta, var_a, var_b, delta, a_gamma, b_gamma, verbose);
  double elbo_init = network.elbo_value;
  network.SetOptimisationPars(tol,n_iter_max,natural_gradient,learning_rate_factor_up,learning_rate_factor_down);
  if (network.verbose) network.Summary();
  network.Optimisation();
  if (network.verbose) network.Summary();
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time,
                             Rcpp::Named("var_alpha_u") = network.var_alpha_u,
                             Rcpp::Named("var_alpha_v") = network.var_alpha_v,
                             Rcpp::Named("var_beta_u") = network.var_beta_u,
                             Rcpp::Named("var_beta_v") = network.var_beta_v,
                             Rcpp::Named("var_lambda") = network.var_lambda,
                             Rcpp::Named("var_delta") = network.var_delta,
                             Rcpp::Named("var_a") = network.var_a,
                             Rcpp::Named("var_b") = network.var_b,
                             Rcpp::Named("learning_rates_u") = network.learning_rates_alpha_beta_u,
                             Rcpp::Named("learning_rates_v") = network.learning_rates_alpha_beta_v,
                             Rcpp::Named("elbo_values") = network.elbo_values_store,
                             Rcpp::Named("elbo_init") = elbo_init,
                             Rcpp::Named("elbo_final") = network.elbo_value));
}
