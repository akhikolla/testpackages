
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

void calc_M1_M2_Hadamard(Eigen::Map<Eigen::ArrayXd> & M1, const Eigen::Map<Eigen::ArrayXd> & M2){ 
 
  //  This function implements:
  //  the Hadamard product $M1 M2$
  //  This function makes in-place computations
  
     M1 *= M2;   
}


