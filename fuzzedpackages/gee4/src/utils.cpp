#include "utils.h"

namespace dragonwell {

  arma::mat VecToUpperTrimatCol(int n, const arma::vec& x, bool diag) {
    arma::mat X = arma::eye<arma::mat>(n, n);

    // make empty matrices
    arma::mat RowIdx(n, n, arma::fill::zeros);
    arma::mat ColIdx(n, n, arma::fill::zeros);

    // fill matrices with integers
    arma::vec idx = arma::linspace<arma::vec>(1, n, n);
    RowIdx.each_col() += idx;
    ColIdx.each_row() += trans(idx);

    // assign upper triangular elements
    // the >= allows inclusion of diagonal elements
    if (diag)
      X.elem(arma::find(RowIdx <= ColIdx)) = x;
    else
      X.elem(arma::find(RowIdx < ColIdx)) = x;

    return X;
  }

  arma::mat VecToLowerTrimatCol(int n, const arma::vec& x, bool diag) {
    arma::mat X = arma::eye<arma::mat>(n, n);

    // make empty matrices
    arma::mat RowIdx(n, n, arma::fill::zeros);
    arma::mat ColIdx(n, n, arma::fill::zeros);

    // fill matrices with integers
    arma::vec idx = arma::linspace<arma::vec>(1, n, n);
    RowIdx.each_col() += idx;
    ColIdx.each_row() += trans(idx);

    // assign upper triangular elements
    // the >= allows inclusion of diagonal elements
    if (diag)
      X.elem(arma::find(RowIdx >= ColIdx)) = x;
    else
      X.elem(arma::find(RowIdx > ColIdx)) = x;

    return X;
  }

  arma::mat ltrimat(arma::uword n, const arma::vec& x, bool diag, bool byrow) {
    arma::mat X;

    if (byrow)
      X = arma::trans(VecToUpperTrimatCol(n, x, diag));
    else
      X = VecToLowerTrimatCol(n, x, diag);

    return X;
  }

  arma::vec UpperTrimatToVecCol(const arma::mat& X, bool diag) {
    int n = X.n_rows;
    arma::vec x;

    // make empty matrices
    arma::mat RowIdx(n, n, arma::fill::zeros);
    arma::mat ColIdx(n, n, arma::fill::zeros);

    // fill matrices with integers
    arma::vec idx = arma::linspace<arma::vec>(1, n, n);
    RowIdx.each_col() += idx;
    ColIdx.each_row() += trans(idx);

    // assign upper triangular elements
    // the >= allows inclusion of diagonal elements
    if (diag)
      x = X.elem(arma::find(RowIdx <= ColIdx));
    else
      x = X.elem(arma::find(RowIdx < ColIdx));

    return x;
  }

  arma::vec LowerTrimatToVecCol(const arma::mat& X, bool diag) {
    int n = X.n_rows;
    arma::vec x;

    // make empty matrices
    arma::mat RowIdx(n, n, arma::fill::zeros);
    arma::mat ColIdx(n, n, arma::fill::zeros);

    // fill matrices with integers
    arma::vec idx = arma::linspace<arma::vec>(1, n, n);
    RowIdx.each_col() += idx;
    ColIdx.each_row() += trans(idx);

    // assign upper triangular elements
    // the >= allows inclusion of diagonal elements
    if (diag)
      x = X.elem(arma::find(RowIdx >= ColIdx));
    else
      x = X.elem(arma::find(RowIdx > ColIdx));

    return x;
  }

  arma::vec lvectorise(const arma::mat& X, bool diag, bool byrow) {
    arma::vec x;

    if (byrow)
      x = UpperTrimatToVecCol(X.t(), diag);
    else
      x = LowerTrimatToVecCol(X, diag);

    return x;
  }

  arma::mat corr_cs(double rho, arma::uword dim) {
    arma::mat I = arma::eye(dim, dim);
    arma::vec one = arma::ones<arma::vec>(dim);
    arma::mat J = one * one.t();
    arma::mat result = (1 - rho) * I + rho * J;
    return result;
  }

  arma::mat corr_ar1(double rho, int dim) {
    arma::mat result = arma::eye(dim, dim);
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        result(i, j) = std::pow(rho, std::abs(i - j));
      }
    }
    return result;
  }

  arma::vec join_vecs(std::initializer_list<arma::vec> vecs) {
    arma::vec result;
    result.reset();

    for (auto iter = vecs.begin(); iter != vecs.end(); ++iter) {
      result.insert_rows(result.n_elem, *iter);
    }

    return result;
  }
}
