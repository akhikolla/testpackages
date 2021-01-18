/*    all.cpp:  1.:Premultiply input samples y by a square root of the Toeplitz  
#               correlation matrix correspdg. to the input partial 
#               autocorrelations partial.
#               2.:Then compute a HAC-type test statistic on each of the input 
#               samples (Equation (7) of Preinerstorfer and Poetscher (2016); 
#               cf. also the discussion in Section 3.2.1 A.); or compute a 
#               weighted Eicker type test statistic on each of the input samples
#               (Eicker 1967, cf. also Section 3.2.1 B of Preinerstorfer and 
#               Poetscher 2016).
#               3.:Return values of test statistic on input samples. 
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
#include "csart.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;

// [[Rcpp::export]]
Eigen::VectorXd testvals(Eigen::Map<Eigen::MatrixXd> y,
                   Eigen::Map<Eigen::MatrixXd> Coefpremult,
                   Eigen::Map<Eigen::VectorXd> partial,
                   Eigen::Map<Eigen::MatrixXd> X,
                   Eigen::Map<Eigen::MatrixXd> Wmat, 
                   Eigen::Map<Eigen::MatrixXd> Bmat,
                   Eigen::Map<Eigen::MatrixXd> R,
                   int dim,
                   int Nrep,
                   int q,
                   int cores, 
                   bool Eicker) {

//Premultiplication step:                
MatrixXd P = premult(partial, dim);
//compute test statistics:

//Preparation of auxiliaries:
MatrixXd obs = P * y;                
MatrixXd Coefs = Coefpremult * obs;                
MatrixXd Rbmat = R * Coefs;
MatrixXd umat = obs - X * Coefs;
VectorXd vals(Nrep);

//Computation step:
if ( Eicker ) {
 vals= ctestE(umat, Rbmat, Wmat, Bmat, cores);
} else {
 vals= ctest(umat, Rbmat, Wmat, Bmat, cores);
}
return vals;
}


