/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _ccaPP_UTILS_H
#define _ccaPP_UTILS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// functions to export to R
RcppExport SEXP R_l1Median(SEXP R_x); // for testing
RcppExport SEXP R_fastMedian(SEXP R_x);
RcppExport SEXP R_fastMAD(SEXP R_x, SEXP R_constant);
RcppExport SEXP R_rank(SEXP R_x);     // for testing

// functions to be used within C++
mat covMCD(const mat& x);
vec l1Median(const mat& x);
double median(const vec& x);
double mad(const vec& x, double& center);
double mad(const vec& x, const double& constant, double& center);
uvec order(const vec& x, const bool& decreasing);
uvec order(const vec& x);
vec rank_ccaPP(const vec& x);

#endif
