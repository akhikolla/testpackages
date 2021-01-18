/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _sparseLTSEigen_FASTLASSO_H
#define _sparseLTSEigen_FASTLASSO_H

#define EIGEN_NO_DEBUG

#include <sparseLTSEigen.h>

using namespace Eigen;

// functions to export to R
RcppExport SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
  	SEXP R_subset, SEXP R_normalize, SEXP R_intercept, SEXP R_eps, 
    SEXP R_useGram);

// functions to be used within C++
void fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
  	const bool& useSubset, const VectorXi& subset, const bool& normalize, 
    const bool& useIntercept, const double& eps, const bool& useGram, 
    const bool& useCrit,
    // intercept, coefficients, residuals and objective function are returned 
    // through the following parameters
    double& intercept, VectorXd& beta, VectorXd& residuals, double& crit);

#endif
