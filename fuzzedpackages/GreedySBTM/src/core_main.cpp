#include <iostream>
#include <RcppArmadillo.h>
#include "core_dsbtm.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_ICLExact(arma::cube adj, arma::mat z, bool verbose)
{
  double computing_time;
  arma::wall_clock timer;
  timer.tic();
  dsbtm network(adj,z,100,verbose);
  computing_time = timer.toc();
  if (verbose) network.Print();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time, 
                             Rcpp::Named("prior_value") = network.prior_value, 
                             Rcpp::Named("likelihood_value") = network.likelihood_value, 
                             Rcpp::Named("icl_value") = network.posterior_value));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_GreedyICL(arma::cube adj, arma::mat z, unsigned int max_n_iter, bool verbose)
{
  double posterior_value_start, posterior_value_end, computing_time;
  arma::wall_clock timer;
  timer.tic();
  dsbtm network(adj,z,max_n_iter,verbose);
  posterior_value_start = network.posterior_value;
  network.GreedyOptimisation();
  posterior_value_end = network.posterior_value;
  computing_time = timer.toc();
  network.DebugCheckAllValues();// only for debug
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time, 
                             Rcpp::Named("icl_start") = posterior_value_start, 
                             Rcpp::Named("icl_trace") = network.greedy_icl_store, 
                             Rcpp::Named("icl_end") = posterior_value_end, 
                             Rcpp::Named("allocations") = network.z));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_GreedyMerge(arma::cube adj, arma::mat z, bool verbose)
{
  double posterior_value_start, posterior_value_end, computing_time;
  arma::wall_clock timer;
  timer.tic();
  dsbtm network(adj,z,100,verbose);
  posterior_value_start = network.posterior_value;
  network.MergeUpdates();
  posterior_value_end = network.posterior_value;
  computing_time = timer.toc();
  return (Rcpp::List::create(Rcpp::Named("computing_time") = computing_time, 
                             Rcpp::Named("icl_start") = posterior_value_start, 
                             Rcpp::Named("icl_end") = posterior_value_end, 
                             Rcpp::Named("allocations") = network.z));
}


