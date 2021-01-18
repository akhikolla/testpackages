#include "choleskyDecomposition.h"

#include <algorithm>

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

double* choleskyDecomposition(const int size, const NumericVector &covariances) {
  const char uplo = 'U';
  const int bands = std::min<int>(covariances.size() - 1, size - 1);
  const int ldA = bands + 1;
  
  double* A;
  A = new double[size * ldA];
  
  for (int i = 0; i <= bands; ++i) {
    for (int j = i; j < size; ++j) {
      A[j * ldA + bands - i] = covariances[i];
    }
  }
  
  int info;  
  F77_CALL(dpbtf2)(&uplo, &size, &bands, A, &ldA, &info);
  
  if (info != 0) {
    if (info < 0) {
      stop("the %d-th argument of the covariance matrix had an illegal value", -info);
    } else {
      stop("a deconvolution could not be performed, since the leading minor of order %d of the covariance matrix is not positive definite. Plese use a different regularization.", info);
    }
  }
  return A;
}
