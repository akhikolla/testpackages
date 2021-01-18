// init.h
#ifndef INIT_H_
#define INIT_H_

#include <vector>
#include <tuple>

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <Rcpp.h>

#include "grid.h"
#include "prep.h"

using prep::Es;

namespace init {
////////////////////////////////////////////////////////////////

// Note: this algorithm is too slow if we use more than two E grids
Rcpp::List gridSearch(Rcpp::NumericMatrix const& exprs,
                      Rcpp::NumericMatrix const& weights,
                      Rcpp::NumericVector const& D,
                      Rcpp::NumericVector const& N8,
                      Rcpp::NumericVector const& Ccos,
                      Rcpp::NumericVector const& Csin,
                      Es const&                  es,
                      std::string&               e) {
  std::size_t const ngenes   = exprs.ncol();
  std::size_t const nsamples = exprs.nrow();
  std::size_t const nparams = 6; // env.e.*, gate.e.*
  std::size_t const ninputs = 9; // 1,D,N8,C.cos,C.sin,D*C.cos,D*C.sin,E(e),D*E(e)
  
  // penalty given by W*(Y-XB)^2
  // - for each gene, w*(y-Xb)^2 wehre w and y are vecs of length nsamples
  Eigen::Map<Eigen::MatrixXd> const Y(REAL(exprs), nsamples, ngenes);
  Eigen::Map<Eigen::MatrixXd> const W(REAL(weights), nsamples, ngenes);
  Eigen::MatrixXd                   X(nsamples, ninputs);
  Eigen::MatrixXd                   B(ninputs, ngenes);
  Eigen::VectorXd                   dev(ngenes); // deviance

  Rcpp::NumericMatrix               best_params(nparams, ngenes);
  Rcpp::NumericVector               best_devs(ngenes, -1); // deviance is never negative
  
  double* x = X.data();
  for (std::size_t i = 0; i < nsamples; ++i, ++x) *x = 1;                       // 1
  for (auto u = D.begin();  u != D.end();  ++u, ++x) *x = *u;                   // D
  for (auto u = N8.begin(); u != N8.end(); ++u, ++x) *x = *u;                   // N8
  for (auto u = Ccos.begin(); u != Ccos.end(); ++u, ++x) *x = *u;               // C.cos
  for (auto u = Csin.begin(); u != Csin.end(); ++u, ++x) *x = *u;               // C.sin
  auto v = D.begin();
  for (auto u = Ccos.begin(); u != Ccos.end(); ++v, ++u, ++x) *x = *v * (*u);   // DC.cos
  v = D.begin();
  for (auto u = Csin.begin(); u != Csin.end(); ++v, ++u, ++x) *x = *v * (*u);   // DC.sin

  double* x0 = x;
  // the parameter grid
  for (auto e0 = es.cbegin(); e0 != es.cend(); ++e0) {
    x = x0;
    for (auto u = e0.cbegin(); u != e0.cend(); ++u, ++x) *x = *u;             // E(e)
    auto v = D.begin();
    for (auto u = e0.cbegin(); u != e0.cend(); ++v, ++u, ++x) *x = *v * (*u); // D*E(e)

    // linear regression with weighted square sum penalties
    for (auto g = 0; g < B.cols(); ++g) {
      Eigen::VectorXd const& w = W.col(g);
      B.col(g) = (X.transpose() * w.asDiagonal() * X).ldlt()
        .solve(X.transpose() * w.cwiseProduct(Y.col(g)));
    }
    
    dev = (Y - X * B).cwiseAbs2().cwiseProduct(W).colwise().sum();
    // B = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
    // dev = (Y - X * B).colwise().squaredNorm();

    double* q = dev.data();
    auto bd = best_devs.begin();
    auto bp = best_params.begin();
    for (std::size_t i = 0; i < ngenes; ++i, ++q, ++bd, bp += nparams){
      if (!Rcpp::NumericVector::is_na(*q) && (*q < *bd || *bd < 0)) {
        *bd = *q;
        auto params = e0.coord();
        // ordering of E coords is dictated by performace.
        // we here reorder them to match the interface with R
        *(bp + 0) = std::get<5>(params); // period.e
        *(bp + 1) = std::get<3>(params); // env.e.amplitude
        *(bp + 2) = std::get<4>(params); // env.e.threshold
        *(bp + 3) = std::get<2>(params); // gate.e.phase
        *(bp + 4) = std::get<0>(params); // gate.e.amplitude
        *(bp + 5) = std::get<1>(params); // gate.e.threshold
      }
    }
  }

  Rcpp::List dimnames = exprs.attr("dimnames");
  // DOC: this order must be consistent with prep::compEs_ (ugly..)
  auto rownames = Rcpp::CharacterVector::create("env."  + e + ".period",
                                                "env."  + e + ".amplitude",
                                                "env."  + e + ".threshold",
                                                "gate." + e + ".phase",
                                                "gate." + e + ".amplitude",
                                                "gate." + e + ".threshold");
  auto colnames = (dimnames.size() == 2) ? dimnames[1] : Rcpp::CharacterVector::create();
  best_params.attr("dimnames") = Rcpp::List::create(rownames, colnames);
  best_devs.attr("names") = colnames;
  return Rcpp::List::create(Rcpp::Named("init.devs") = best_devs,
                            Rcpp::Named("init.params") = best_params);
}

////////////////////////////////////////////////////////////////
} // namespace init
#endif // INIT_H_
