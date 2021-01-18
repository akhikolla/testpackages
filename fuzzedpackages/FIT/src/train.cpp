// train.cpp
// library(Rcpp)
// Sys.setenv(PKG_CXXFLAGS = "-Wall -g -O3 -std=c++11 -I/Users/yuri/.R/RcppEigen/include")
// sourceCpp('./lib/FIT/src/train.cpp')
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Rcpp.h>
#include <RcppEigen.h>

#include "grid.h"
#include "model.h"
#include "prep.h"
#include "init.h"
#include "optim.h"

using prep::Es;

// [[Rcpp::depends(RcppEigen)]]

// XXX: at the moment, assumes that envs is a list of singletons

////////////////////////////////////////////////////////////////
// DOC
//
// Interface: Communicates with the calling R process by names
// - Following names are a part of our protocol:
//   - gene names
//   - attribute$data names: "times.pickup", "times.of.day", "age", "type"
//   - weather$data names: "times.of.day" (consistent with that in attribute$data)
//                         "wind", "temperature", "humidity", "atmosphere",
//                         "precipitation", "radiation"
//   - grid coordinate names: for e in weather_factor_names
//                              env.e.period, env.e.amplitude, env.e.threshold,
//                              gate.e.phase, gate.e.amplitude, gate.e.threshold
// - Changing these names in attribute$data etc. *may break* the C++ code
//   (although at the moment we only use the coord names
//    and the weather factor names (the latter only implicitly))
//
// Args:
// - exprs: Expression data packed in an R matrix of dim [samples.n * genes.n]
// - attribute_data: An R data frame with samples.n obs (cols: times.of.day etc.)
// - weather_data:   An R data frame with long enough obs to cover the periods
//                   of interest as are computed from attribute$data$times.pickup
//                   and grid_coordinates$gate.e.period
// - env_combinations: A list of envs (CharacterVectors) where each vector holds
//                     the names of weather factors (see above) that are to
//                     be included in an env
// - grid_coordinates: A list of coordinates (NumericVectors or IntegerVectors)

// gridSearch

// [[Rcpp::export]]
Rcpp::List initParamsAndDevs(Rcpp::NumericMatrix const   exprs,
                             Rcpp::NumericMatrix const   weights,
                             Rcpp::DataFrame const       attribute_data,
                             Rcpp::DataFrame const       weather_data,
                             Rcpp::CharacterVector const env_factors,
                             Rcpp::List const            grid_coordinates,
                             Rcpp::IntegerVector const   data_step,
                             Rcpp::IntegerVector const   time_step) {
  // XXX:change env_factors to a list of weather_factor -> genes
  // XXX:revise this invariant:
  if (exprs.nrow() != attribute_data.nrows())
    throw Rcpp::exception("nrows of expr and attribute_data do not match.");
  if (data_step.size() != 1 || time_step.size() != 1)
    throw Rcpp::exception("data_step and time_step are supposed be scalars.");

  Rcpp::Rcout << "# Prep (grids)\n";
  Rcpp::Rcout << "# - D, type, C\n";
  Rcpp::NumericVector const D  = attribute_data["age"];
  Rcpp::NumericVector const N8 = attribute_data["type"];

  Rcpp::NumericVector timesOfDay = attribute_data["times.of.day"];
  Rcpp::NumericVector Ccos = model::clockCos(timesOfDay);
  Rcpp::NumericVector Csin = model::clockSin(timesOfDay);

  std::vector<std::unique_ptr<Es>> es(env_factors.size());
  auto p = es.begin();
  
  for (auto& e_ : env_factors) {
    auto e = Rcpp::as<std::string>(e_);
    Rcpp::Rcout << "# - E(" + e + ")\n";
    *p++ = prep::makeEs(attribute_data["times.pickup"],
                        attribute_data["times.of.day"],
                        weather_data[e],
                        // XXX: as IntegerVector?
                        grid_coordinates["env." + e + ".period"],
                        grid_coordinates["env." + e + ".amplitude"],
                        grid_coordinates["env." + e + ".threshold"],
                        grid_coordinates["gate." + e + ".phase"],
                        grid_coordinates["gate." + e + ".amplitude"],
                        grid_coordinates["gate." + e + ".threshold"],
                        data_step[0], time_step[0]);
  }
  
  Rcpp::Rcout << "# Init (grid search)\n";
  // fix params for e from (D,N8,C,DC,E(e),DE(e)) for a given e
  // XXX: How to fix params for C when an env contains more than one e
  // - The simplest heuristic is to use the value for the E that gave the best deviance.
  Rcpp::List res(env_factors.size());
  res.attr("names") = env_factors;
  auto u = res.begin();
  auto v = es.cbegin();
  for (auto& e_ : env_factors) { // XXX: actually should be envs here
    auto e = Rcpp::as<std::string>(e_);
    Rcpp::Rcout << "# - init params for " + e << '\n';
    // XXX: make matrix exprs
    // XXX: update initParamsAndDev1()
    *u = init::gridSearch(exprs, weights, D, N8, Ccos, Csin, *(*v), e);
    ++u; ++v;
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix inputVars(Rcpp::NumericVector const   params,
                              Rcpp::CharacterVector const env,
                              Rcpp::DataFrame const       attribute_data,
                              Rcpp::DataFrame const       weather_data,
                              Rcpp::IntegerVector const   data_step,
                              Rcpp::IntegerVector const   time_step) {
  if (data_step.size() != 1 || time_step.size() != 1)
    throw Rcpp::exception("data_step and time_step are supposed be scalars.");
  std::size_t const nsamples = attribute_data.nrows();
  std::size_t const ninputs = 7 + 2 * env.size(); // 1,D,N8,C.cos,C.sin,DC.cos,DC.sin,E(e),DE(e)
  
  std::unique_ptr<Eigen::MatrixXd> p = optim::inputVars(nsamples, ninputs, params, env,
                                                        attribute_data, weather_data,
                                                        data_step[0], time_step[0]);
  return Rcpp::wrap(*p);
}

// [[Rcpp::export]]
Rcpp::NumericVector devLm(Rcpp::NumericVector const   params,
                          Rcpp::CharacterVector const env,
                          Rcpp::NumericVector const   expr,
                          Rcpp::NumericVector const   weight,
                          Rcpp::DataFrame const       attribute_data,
                          Rcpp::DataFrame const       weather_data,
                          Rcpp::IntegerVector const   data_step,
                          Rcpp::IntegerVector const   time_step) {
  if (data_step.size() != 1 || time_step.size() != 1)
    throw Rcpp::exception("data_step and time_step are supposed be scalars.");
  
  std::size_t const nsamples = expr.size();
  std::size_t const ninputs = 7 + 2 * env.size(); // 1,D,N8,C.cos,C.sin,DC.cos,DC.sin,E(e),DE(e)
  
  std::unique_ptr<Eigen::MatrixXd> p = optim::inputVars(nsamples, ninputs, params, env,
                                                        attribute_data, weather_data,
                                                        data_step[0], time_step[0]);
  // y ~ X.b
  Eigen::Map<Eigen::VectorXd> const y(REAL(expr), nsamples);
  Eigen::MatrixXd const&            X = *p;
  Eigen::VectorXd                   b(ninputs);
  Eigen::Map<Eigen::VectorXd> const w(REAL(weight), nsamples);

  b = (X.transpose() * w.asDiagonal() * X).ldlt()
    .solve(X.transpose() * w.cwiseProduct(y));
  
  double dev = (y - X * b).cwiseAbs2().dot(w);
  // b = (X.transpose() * X).ldlt().solve(X.transpose() * y);
  // double dev = (y - X * b).squaredNorm();
  
  return Rcpp::wrap(dev);
}

// [[Rcpp::export]]
Rcpp::NumericVector coefsLm(Rcpp::NumericVector const   params,
                            Rcpp::CharacterVector const env,
                            Rcpp::NumericVector const   expr,
                            Rcpp::NumericVector const   weight,
                            Rcpp::DataFrame const       attribute_data,
                            Rcpp::DataFrame const       weather_data,
                            Rcpp::IntegerVector const   data_step,
                            Rcpp::IntegerVector const   time_step) {
  if (data_step.size() != 1 || time_step.size() != 1)
    throw Rcpp::exception("data_step and time_step are supposed be scalars.");
  std::size_t const nsamples = expr.size();
  std::size_t const ninputs = 7 + 2 * env.size(); // 1,D,N8,C.cos,C.sin,DC.cos,DC.sin,E(e),DE(e)

  std::unique_ptr<Eigen::MatrixXd> p = optim::inputVars(nsamples, ninputs, params, env,
                                                        attribute_data, weather_data,
                                                        data_step[0], time_step[0]);
  // minimize w*(y-Xb)^2
  // b = (Xt . wX)^{-1}(Xt.wy)
  Eigen::Map<Eigen::VectorXd> const y(REAL(expr), nsamples);
  Eigen::Map<Eigen::VectorXd> const w(REAL(weight), nsamples);
  Eigen::MatrixXd const&            X = *p;
  Eigen::VectorXd                   b(ninputs);

  b = (X.transpose() * w.asDiagonal() * X).ldlt()
    .solve(X.transpose() * w.cwiseProduct(y));
  // b = (X.transpose() * X).ldlt().solve(X.transpose() * y);

  return Rcpp::wrap(b);
}

// Walk around Mingw g++/Rcpp::export bug
// - the name of the last alphabetical entry with RcppEigen gets mungled
//   in a strange way on Windows.
// [[Rcpp::export]]
Rcpp::NumericMatrix zzzRcppExportBug() { Eigen::MatrixXd X(10,10); return Rcpp::wrap(X); }
