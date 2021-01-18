
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

  Eigen::MatrixXd calc_M1timesM2v( const Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::MatrixXd> & M2, const Eigen::Map<Eigen::ArrayXd> & v){

  //  This function implements:
  //  M1 %*% (t(M2) * v) 

  const unsigned int m = M2.rows();
  Eigen::MatrixXd Intermediate =  M2;
  for (unsigned int i = 0; i < m; ++i){
    Intermediate.row(i) = M2.row(i).transpose().array() * v;  
  }  
  
  return ( M1 * Intermediate.transpose() );  
}

 


