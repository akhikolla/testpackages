/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _sparseLTSEigen_UTILS_H
#define _sparseLTSEigen_UTILS_H

#define EIGEN_NO_DEBUG

#include <sparseLTSEigen.h>

using namespace Eigen;

// functions to export to R
RcppExport SEXP R_findSmallest(SEXP R_x, SEXP R_h);
RcppExport SEXP R_partialOrder(SEXP R_x, SEXP R_h);
//SEXP R_partialSort(SEXP R_x, SEXP R_h);

// functions to be used within C++
VectorXi findSmallest(const VectorXd& x, const int& h);
VectorXi partialOrder(const VectorXd& x, const int& h);
//VectorXi order(const VectorXd& x);

#endif
