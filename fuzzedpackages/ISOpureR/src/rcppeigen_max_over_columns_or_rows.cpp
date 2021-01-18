// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_max_over_columns_or_rows(const Eigen::Map<Eigen::MatrixXd> & A, int & dimen){
	Eigen::MatrixXd Amax;
	if (dimen == 1){
		Amax = A.colwise().maxCoeff();
	}
	else if (dimen ==2){
		Amax = A.rowwise().maxCoeff();
	}
	return Amax;
}
