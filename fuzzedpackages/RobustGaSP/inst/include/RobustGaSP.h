/*
// RobustGaSP.h
// Jesus Palomo and and Mengyang Gu, February 2016
//
// Copyright (C)  2016 Jesus Palomo and Mengyang Gu
//
// This file is part of the Robust GaSP Package available at CRAN.
//
// Robust GaSP is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Robust GaSP is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/

#ifndef _robustgasp_ctools_H
#define _robustgasp_ctools_H
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

typedef   Eigen::VectorXi        iVec;
typedef   Eigen::Map<iVec>      MapiVec;
typedef   Eigen::MatrixXd         Mat;
typedef   Eigen::Map<Mat>        MapMat;
typedef   Eigen::VectorXd         Vec;
typedef   Eigen::Map<Vec>        MapVec;
typedef   Eigen::ArrayXd          Ar1;
typedef   Eigen::Map<Ar1>        MapAr1;
typedef   Eigen::ArrayXXd         Ar2;
typedef   Eigen::Map<Ar2>        MapAr2;

using      Rcpp::List;
using      Rcpp::Named;
using      Rcpp::NumericVector;
using      Rcpp::as;
using      Rcpp::wrap;

typedef std::vector<std::map<std::string, int> > mp_container;


NumericMatrix matern_2_5 (ListMatrix R0, NumericVector beta);
//inline static Mat matern_5_2_funct (const MapMat &d, double beta_i);
Mat matern_2_5_funct (const MapMat & d, double beta_i);
RcppExport SEXP check_openmp();

#endif
