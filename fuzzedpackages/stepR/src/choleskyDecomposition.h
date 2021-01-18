#ifndef STEPR_H_CHOL
#define STEPR_H_CHOL

#include <Rcpp.h>

using namespace Rcpp;

// cholesky decomposition, returned variable has to be deleted in the later code
double* choleskyDecomposition(const int size, const NumericVector &covariances);

#endif
