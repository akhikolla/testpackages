/*    test.statistic.cpp:  Used in the computation of a HAC-type test statistic 
#               (Equation (7) of Preinerstorfer and Poetscher (2016); cf. also 
#               the discussion in Section 3.2.1 A.) in function testvals
#               defined in all.cpp.
#               
#    Copyright (C) 2016  David Preinerstorfer
#    david.preinerstorfer@econ.au.dk
#
#    This file is a part of acrt.
#
#    acrt is free software; you can redistribute it and/or modify
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

#ifdef _OPENMP
# include <omp.h>
#endif
#include <RcppEigen.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

// [[Rcpp::export]]

Eigen::VectorXd ctest(Eigen::MatrixXd umat,
                Eigen::MatrixXd Rbmat,
                Eigen::MatrixXd Wmat, 
                Eigen::MatrixXd Bmat,
                int cores) {

int Nrep = umat.cols(), q = Bmat.rows();
VectorXd testevals(Nrep);

#ifdef _OPENMP
omp_set_num_threads(cores);
#pragma omp parallel for
#endif
for(int i = 0; i < Nrep; i++){
MatrixXd Wmat0 = Bmat * umat.col(i).asDiagonal();
MatrixXd Wmat1(q, q);
 for(int k = 0; k < q; k++){
 VectorXd v =  Wmat * Wmat0.row(k).transpose();
  for(int l = 0; l <=k; l++){
  Wmat1(k, l) =  Wmat0.row(l) * v;
  }
 }
FullPivLU<MatrixXd>luWmat1(Wmat1.selfadjointView<Lower>());
if(luWmat1.isInvertible()){
testevals[i] = Rbmat.col(i).transpose() * luWmat1.solve(Rbmat.col(i));
} else {
testevals[i] = 0;
}
}
return testevals;
}
