// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include "hera/wasserstein/include/wasserstein.h"


// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double wasserstein_distance(const Rcpp::NumericMatrix & Diag1, const Rcpp::NumericMatrix & Diag2, int q,  double internal_p, double delta = 0.01, double initial_eps = 1, double eps_factor=5.0)
{
  hera::AuctionParams<double> params;
  params.wasserstein_power = q;
  params.delta = delta;
  params.internal_p = internal_p;
  if (initial_eps != 0)
    params.initial_epsilon = initial_eps ;
  
  if (eps_factor != 0.)
    params.epsilon_common_ratio = eps_factor;
  
  using PairVector = std::vector<std::pair<double, double>>;
  PairVector diagramA, diagramB;
  
// itera su tutto A e tutto B
  for (int i = 0; i < Diag1.nrow(); i++)
    diagramA.push_back(std::make_pair(Diag1(i,0),Diag1(i,1) ));
  
  for (int j = 0; j < Diag2.nrow(); j++)
    diagramB.push_back(std::make_pair(Diag2(j,0),Diag2(j,1) ));

  return hera::wasserstein_dist(diagramA, diagramB, params);
}

