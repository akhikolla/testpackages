/*    premult.cpp: Compute a (lower triangular) square root of the dim x dim - 
#               dimensional AR(p) correlation matrix correspdg. to the input  
#               partial autocorrelation sequence partial.
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

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

/* Construct Premultiplication Matrix Premult: Premult * t(Premult) is the 
dim x dim - dimensional Toeplitz correlation matrix with partial 
autocorrelations partial */

Eigen::MatrixXd premult(Eigen::VectorXd partial, int dim){

int p = partial.size();
VectorXd d(p+1);
VectorXd dvec = VectorXd::Constant(dim, 1); 
MatrixXd alpha = MatrixXd::Constant(p+1, p+1, 0); 
MatrixXd Premult = MatrixXd::Constant(dim, dim, 0);
d(0) = 1; 
alpha(0, 0) = 1;

/*Durbin-Levinson recursion to obtain prediction coefficients from partial
autocorrelations:*/

for(int n = 0; n < p; n++){
 d(n+1) = d(n)*(1-std::pow(partial(n), 2));
 alpha(n+1, 0) = - partial(n);
 alpha(n+1, n+1) = 1;
 if( n >= 1 ){
  for(int j = 1; j <= n; j++){
  alpha(n+1, j) = alpha(n, j-1) - partial(n)*alpha(n, n-j);
  }
 }
}

/*Construction of the square root of the correlation matrix:*/

if(dim >= p+1){
 Premult.topLeftCorner(p+1, p+1) = alpha;
 for(int i = p+1; i < dim; i++){
  Premult.block(i, i-p, 1, p+1) = alpha.row(p);
 }
 dvec = std::sqrt(d(p)) * dvec;
 dvec.head(p+1) = d.array().sqrt(); 
} else {
 Premult = alpha.topLeftCorner(dim, dim);
 dvec = d.head(dim).array().sqrt();
}
MatrixXd Pmat = Premult.triangularView<Lower>().solve(MatrixXd::Identity(dim, dim)) * dvec.asDiagonal(); 
return Pmat;                                      
}
