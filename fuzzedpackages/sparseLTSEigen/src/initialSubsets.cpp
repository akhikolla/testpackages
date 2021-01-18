/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include "initialSubsets.h"

using namespace Rcpp;
using namespace Eigen;

// function used internally, which computes lasso fits for subsets containing a 
// small number of observations (typically only 3) and returns the indices of 
// the respective h observations with the smallest absolute residuals
MatrixXi sparseSubsets(const MatrixXd& x, const VectorXd& y,
  	const double& lambda, const int& h, const MatrixXi& subsets, 
    const bool& normalize, const bool& useIntercept, 
    const double& eps, const bool& useGram) {
	const int nsamp = subsets.cols();
	MatrixXi indices(h, nsamp);
	for(int k = 0; k < nsamp; k++) {
		// compute lasso fit
  	double intercept, crit;
		VectorXd coefficients, residuals;
    fastLasso(x, y, lambda, true, subsets.col(k), normalize, useIntercept, 
        eps, useGram, false, intercept, coefficients, residuals, crit);
		// find h observations with smallest absolute residuals
		indices.col(k) = findSmallest(residuals.cwiseAbs(), h);
	}
	return indices;
}

// R interface to sparseSubsets()
SEXP R_sparseSubsets(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_h, 
    SEXP R_subsets, SEXP R_normalize, SEXP R_intercept, SEXP R_eps, 
    SEXP R_useGram) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);              // predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);  // reuse memory
	NumericVector Rcpp_y(R_y);              // response
	Map<VectorXd> y(Rcpp_y.begin(), n);     // reuse memory
  double lambda = as<double>(R_lambda);
  int h = as<int>(R_h);
	IntegerMatrix Rcpp_subsets(R_subsets);  // subset to use for computation
  const int s = Rcpp_subsets.nrow(), nsamp = Rcpp_subsets.ncol();
	MatrixXi subsets(s, nsamp);
	for(int j = 0; j < nsamp; j++) {
		for(int i = 0; i < s; i++) {
			subsets(i,j) = Rcpp_subsets(i,j) - 1;
		}
	}
  bool normalize = as<bool>(R_normalize);
	bool useIntercept = as<bool>(R_intercept);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results
	MatrixXi indices = sparseSubsets(x, y, lambda, h, subsets,
			normalize, useIntercept, eps, useGram);
	IntegerVector Rcpp_indices = wrap(indices);
  Rcpp_indices = Rcpp_indices + 1;
  return Rcpp_indices;
}
