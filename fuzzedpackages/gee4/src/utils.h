#ifndef UTILS_H_
#define UTILS_H_

#include <RcppArmadillo.h>

namespace dragonwell {

  // vector to upper triangular matrix with column-major order (Fortran style)
  // arma::mat VecToUpperTrimatCol(int n, const arma::vec& x, bool diag = false);

  // vector to lower triangular matrix with column-major order (Fortran style)
  // arma::mat VecToLowerTrimatCol(int n, const arma::vec& x, bool diag = false);

  arma::mat ltrimat(arma::uword n, const arma::vec& x, bool diag = false,
                    bool byrow = true);

  // upper triangular matrix to vector with column-major order (Fortran style)
  // arma::vec UpperTrimatToVecCol(const arma::mat& X, bool diag = false);

  // lower triangular matrix to vector with column-major order (Fortran style)
  // arma::vec LowerTrimatToVecCol(const arma::mat& X, bool diag = false);

  arma::vec lvectorise(const arma::mat& X, bool diag = false, bool byrow = true);

  arma::mat corr_cs(double rho, arma::uword dim);
  arma::mat corr_ar1(double rho, int dim);

  arma::vec join_vecs(std::initializer_list<arma::vec> vecs);
}

#endif  // UTILS_H_
