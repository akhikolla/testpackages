//  arma_util.cpp: implementation of some utility functions for extension of
//                 Armadillo library
//  This file is part of jmcm.
//
//  Copyright (C) 2015-2016 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/

#include "arma_util.h"

namespace pan {

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

arma::mat ltrimat(int n, const arma::vec& x, bool diag, bool byrow) {
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
}  // namespace pan
