/*    test.statistic.cpp:  
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
                Eigen::VectorXd Wvec, 
                Eigen::MatrixXd Bmat,
                int cores,
                double tol) {

int Nrep = umat.cols(), q = Bmat.rows();
VectorXd testevals(Nrep), sqrtwvec = Wvec.cwiseSqrt();


if(q > 1){ 

  MatrixXd Bmatwsqrt = Bmat * sqrtwvec.asDiagonal();
  #ifdef _OPENMP
  omp_set_num_threads(cores);
  #pragma omp parallel for
  #endif
  for(int i = 0; i < Nrep; i++){
  MatrixXd Wmat0 = Bmatwsqrt * umat.col(i).asDiagonal();
  MatrixXd Wmat1(q, q); 
  Wmat1.setZero().selfadjointView<Lower>().rankUpdate(Wmat0);
  FullPivLU<MatrixXd>luWmat1(Wmat1.selfadjointView<Lower>());
  if(tol > 0){
  luWmat1.setThreshold(tol);
  }
  if(luWmat1.isInvertible()){
  testevals[i] = Rbmat.col(i).transpose() * luWmat1.solve(Rbmat.col(i));
  } else {
  testevals[i] = -1;
  }
  }

} else {

  VectorXd Bmatwsqrt = sqrtwvec.cwiseProduct(Bmat.transpose()); 
  MatrixXd Wmat = Bmatwsqrt.asDiagonal() * umat;
  testevals = Wmat.colwise().squaredNorm(); 	
  VectorXd Rbsquare = Rbmat.colwise().squaredNorm();
  if(tol <= 0){
  tol = DBL_EPSILON;
  }
  #ifdef _OPENMP
  omp_set_num_threads(cores);
  #pragma omp parallel for
  #endif
  for(int i = 0; i < Nrep; i++){
  if(testevals[i] < tol){
  testevals[i] = -1;
  } else {
  testevals[i] = Rbsquare[i]/testevals[i];
  }
  }
} 

return testevals;
}
