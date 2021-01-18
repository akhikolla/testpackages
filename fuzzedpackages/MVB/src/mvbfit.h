//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include "lps.h"
#include "comb.h"
#include "gme.h"

// declaration of functions

// collect resutls for fit and lps
void inline collectFIT(Rcpp::List&, lps::lps&);
void inline collectLPS(Rcpp::List&, lps::lps&);
// tune method string to integer
int getMethod(const std::string&);
// calculate eS from fx
RcppExport SEXP get_eS(SEXP, SEXP);
// evalate logarithm of likelihood of the model
RcppExport SEXP loglike(SEXP, SEXP, SEXP, SEXP);
// fit multivariate Bernoulli model
RcppExport SEXP mvbfit(SEXP, SEXP, SEXP, SEXP, SEXP);
// fit multivariate Bernoulli with l1 penalty
RcppExport SEXP mvblps(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// fit multivairate Bernoulli mixed effect model
RcppExport SEXP mvbme(SEXP,  SEXP, SEXP, SEXP, SEXP, SEXP);
// fit univariate model
RcppExport SEXP unifit(SEXP, SEXP, SEXP, SEXP);
// fit lasso pattern search for univariate
RcppExport SEXP unilps(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// fit multivariate Bernoulli with heuristic order selection
RcppExport SEXP stepfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);



