// optim.h
#ifndef OPTIM_H_
#define OPTIM_H_

// #include <iostream>
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include "grid.h"
#include "prep.h"

using prep::Es;
namespace optim {
////////////////////////////////////////////////////////////////
// XXX: wrap by try to catch index out of bounds? params[] can throw
// ninputs == 5 + 2 * env.size(); // 1,D,N8,C,DC,E(e),DE(e)
std::unique_ptr<Eigen::MatrixXd> inputVars(std::size_t                 nsamples,
                                           std::size_t                 ninputs,
                                           Rcpp::NumericVector const   params,
                                           Rcpp::CharacterVector const env,
                                           Rcpp::DataFrame const       attribute_data,
                                           Rcpp::DataFrame const       weather_data,
                                           int dataStep = 1, int timeStep = 1) {
  if (nsamples != attribute_data.nrows())
    throw Rcpp::exception("nsamples and attribute_data do not match.");
  // env.e.{period,amp,th}, gate.e.{phase,amp,th}
  if (params.size() != (3+3) * env.size())
    throw Rcpp::exception("params.size() and env.size() are inconsistent.");
  
  Rcpp::NumericVector const D  = attribute_data["age"];
  Rcpp::NumericVector const N8 = attribute_data["type"];

  Rcpp::NumericVector timesOfDay = attribute_data["times.of.day"];
  Rcpp::NumericVector Ccos = model::clockCos(timesOfDay);
  Rcpp::NumericVector Csin = model::clockSin(timesOfDay);

  std::vector<std::unique_ptr<Es>> E(env.size());
  auto p = E.begin();
  for (auto const& e_ : env) {
    auto e = Rcpp::as<std::string>(e_);
    *p++ = prep::makeE(attribute_data["times.pickup"],
                       attribute_data["times.of.day"],
                       weather_data[e],
                       std::round(params["env." + e + ".period"]),
                       params["env." + e + ".amplitude"],
                       params["env." + e + ".threshold"],
                       std::round(params["gate." + e + ".phase"]),
                       params["gate." + e + ".amplitude"],
                       params["gate." + e + ".threshold"],
                       dataStep, timeStep);
  }

  std::unique_ptr<Eigen::MatrixXd> X { new Eigen::MatrixXd(nsamples, ninputs) };
  double* x = X->data();
  for (std::size_t i = 0; i < nsamples; ++i, ++x) *x = 1;                     // 1
  for (auto u = D.begin();  u != D.end();  ++u, ++x) *x = *u;                 // D
  for (auto u = N8.begin(); u != N8.end(); ++u, ++x) *x = *u;                 // N8

  for (auto u = Ccos.begin(); u != Ccos.end(); ++u, ++x) *x = *u;             // Ccos
  for (auto u = Csin.begin(); u != Csin.end(); ++u, ++x) *x = *u;             // Csin
  auto v = D.begin();
  for (auto u = Ccos.begin(); u != Ccos.end(); ++u, ++x) *x = *v * (*u);      // D*Ccos
  v = D.begin();
  for (auto u = Csin.begin(); u != Csin.end(); ++u, ++x) *x = *v * (*u);      // D*Csin

  for (auto p = E.cbegin(); p != E.cend(); ++p) {
    auto e0 = (*p)->cbegin(); // E(e) grid is just a point
    for (auto u = e0.cbegin(); u != e0.cend(); ++u, ++x) *x = *u;             // E(e)
    auto v = D.begin();
    for (auto u = e0.cbegin(); u != e0.cend(); ++v, ++u, ++x) *x = *v * (*u); // D*E(e)
  }
  return X;
}

////////////////////////////////////////////////////////////////
} // namespace optim
#endif // OPTIM_H_
