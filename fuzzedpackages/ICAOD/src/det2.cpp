#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

double det2(const Eigen::Map<Eigen::MatrixXd> mat , const bool logarithm = false){
  double det_FIM;

   det_FIM= mat.determinant();

   if(logarithm == true){
    if(det_FIM > 0){
      det_FIM = log(det_FIM);}else
      det_FIM = R_NegInf;
    }
  return det_FIM;
}

