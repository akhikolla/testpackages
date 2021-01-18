
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]] 
 
Eigen::ArrayXd calc_tapply_vect_sum(const Eigen::Map<Eigen::ArrayXd> & v1, const Eigen::Map<Eigen::ArrayXi> & v2){

  //  This function implements:
  //  the equivalent of tapply(v_i1, v_i2, sum) for vectors
   
  unsigned int N = v2.size();
  unsigned int a = v2.maxCoeff();
  Eigen::ArrayXd v3 = Eigen::ArrayXd::Zero(a+1);  
  
  for( unsigned int i=0; i !=N; i++){  
    v3(v2(i)) += v1(i);    
  } 

  return( v3 );  
}
