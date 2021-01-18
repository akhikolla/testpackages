#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::MatrixXd calc_VY( const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::MatrixXd> & A, const double b){ 

  //  This function implements:
  //  as.matrix(M %*% A %*% t(M) + b * diag(1,ni[i]))

  Eigen::MatrixXd M1 = M*A*M.transpose();
  M1.diagonal().array() += b;
  return M1;
}


