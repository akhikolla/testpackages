
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]


void calc_v_a(Eigen::Map<Eigen::ArrayXd> & v, const double & a){

  //  This function implements:
  //  the vector multiplication $a v$ a being a scalar
  //  This function makes in-place computations 

  v *= a;  

}
