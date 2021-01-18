/*    rrev.cpp: 
#
#    Copyright (C) 2020  David Preinerstorfer
#    david.preinerstorfer@ulb.ac.be
#
#    This file is a part of wbsd.
#
#    wbsd is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/
*/

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

// [[Rcpp::export]]

int rrank(Eigen::MatrixXd A, double tol) {
FullPivLU<MatrixXd> lu(A);

if(tol > 0){
 lu.setThreshold(tol);
}

return lu.rank();
}

// [[Rcpp::export]]

Eigen::MatrixXd rkernel(Eigen::MatrixXd A, double tol) {
FullPivLU<MatrixXd> lu(A);

if(tol > 0){
 lu.setThreshold(tol);
}

return lu.kernel();
}
