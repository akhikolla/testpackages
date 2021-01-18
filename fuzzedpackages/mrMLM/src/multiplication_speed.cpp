// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::export]]
SEXP multiplication_speed(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
