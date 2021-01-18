
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Eigen::VectorXd calc_M_v( const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::MatrixXd> & M){ 

  //  This function implements:
  //  M v

  return(M * v);  
}


