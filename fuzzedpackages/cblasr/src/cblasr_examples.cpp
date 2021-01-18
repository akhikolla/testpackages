#include <Rcpp.h>
#include "cblas.h"

//'@title Example: Matrix multiplication using cblas_dgemm
//'@description Matrix multiplication using cblas_dgemm
//'@examples
//'## Expected output:
//'## [ 367.76, 368.12
//'##   674.06, 674.72 ]
//'example_cblas_dgemm()
//'@export
// [[Rcpp::export]]
void example_cblas_dgemm() {
  int lda = 3;

  double A[] = {0.11, 0.12, 0.13, 0.21, 0.22, 0.23};

  int ldb = 2;

  double B[] = {1011, 1012, 1021, 1022, 1031, 1032};

  int ldc = 2;

  double C[] = {0.00, 0.00, 0.00, 0.00};

  /* Compute C = A B */

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 3, 1.0, A, lda,
              B, ldb, 0.0, C, ldc);

  Rcpp::Rcout << "[ " << C[0] << ", " << C[1] << std::endl;
  Rcpp::Rcout << "  " << C[2] << ", " << C[3] << " ]" << std::endl;
}
