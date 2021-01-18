//  arma_util.h: Some utility functions for extension of Armadillo library
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

#ifndef JMCM_SRC_ARMA_UTIL_H_
#define JMCM_SRC_ARMA_UTIL_H_

#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

namespace pan {

// vector to upper triangular matrix with column-major order (Fortran style)
arma::mat VecToUpperTrimatCol(int n, const arma::vec& x, bool diag = false);

// vector to lower triangular matrix with column-major order (Fortran style)
arma::mat VecToLowerTrimatCol(int n, const arma::vec& x, bool diag = false);

arma::mat ltrimat(int n, const arma::vec& x, bool diag = false,
                  bool byrow = true);

// upper triangular matrix to vector with column-major order (Fortran style)
arma::vec UpperTrimatToVecCol(const arma::mat& X, bool diag = false);

// lower triangular matrix to vector with column-major order (Fortran style)
arma::vec LowerTrimatToVecCol(const arma::mat& X, bool diag = false);

arma::vec lvectorise(const arma::mat& X, bool diag = false, bool byrow = true);
}  // namespace pan

#endif  // JMCM_SRC_ARMA_UTIL_H_
