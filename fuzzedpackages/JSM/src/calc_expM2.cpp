
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
void calc_expM2(Eigen::Map<Eigen::ArrayXd> & A){ 

  //  This function implements:
  //  A = exp(A);
  //  This function makes in-place computations

  A = (A.exp()); 
}

