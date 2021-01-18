#ifndef RUNTEST_H_
#define RUNTEST_H_

#include <RcppEigen.h>
#include <vector>
#include <set>
#include "test.h"

std::vector<std::vector<double> > convertToVVd(SEXP const& rMatrix);
std::vector<std::vector<int> > convertToVVi(SEXP const& rMatrix);
std::vector<Rank> downUniVariateRank(Rcpp::NumericMatrix XR);

RcppExport SEXP kullback(SEXP m,SEXP  mu1,SEXP mu2,SEXP  p1, SEXP  p2,SEXP proportion1,SEXP proportion2);
RcppExport SEXP adkhi2(SEXP donnees, SEXP p, SEXP proportion, SEXP mu, SEXP nBootstrap);
RcppExport SEXP adkhi2partial(SEXP donnees, SEXP p, SEXP proportion, SEXP mu, SEXP nBootstrap);
#endif /* RUNTEST_H_ */

