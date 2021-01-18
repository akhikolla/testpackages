#ifndef RUNFUNCTIONS_H_
#define RUNFUNCTIONS_H_

#include <RcppEigen.h>
#include <vector>
#include <set>


std::vector<std::vector<std::vector<int> > > numMat2vvvInt(Rcpp::NumericMatrix XR,std::vector<int> const& m);
RcppExport SEXP freqMultiR(SEXP X,SEXP m);
RcppExport SEXP simulISRR(SEXP n,SEXP m,SEXP mu,SEXP p);
RcppExport SEXP loglikelihood(SEXP X,SEXP mu,SEXP p, SEXP proportion,SEXP m, SEXP iterL, SEXP burnL, SEXP IC, SEXP nb_cpus);
RcppExport SEXP computeProba(SEXP X,SEXP mu,SEXP pi,SEXP m);

#endif /* RUNFUNCTIONS_H_ */

