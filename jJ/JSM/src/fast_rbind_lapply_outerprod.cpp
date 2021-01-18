
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd fast_rbind_lapply_outerprod(Rcpp::List const input){ 

  //  This function implements:
  //  do.call(rbind, lapply(1:n, function(i) apply(t(bi.st[[i]]), 1, function(x) x %o% x)))
          
const unsigned int l = input.size();

Rcpp::NumericMatrix xx = input[1];

const unsigned int nr = xx.nrow();
const unsigned int nc = xx.ncol();

Eigen::MatrixXd U(nr*nr*l, nc);
Eigen::MatrixXd u1 =  Eigen::MatrixXd::Zero( nr, nc );  
Eigen::MatrixXd u2 =  Eigen::MatrixXd::Zero( nr, nr );  

for (unsigned int i = 0; i != l; ++i){  
  u1 =  input[i];
  for (unsigned int j = 0; j != nc; ++j){
    u2 = u1.col(j) *  u1.col(j).adjoint();
    U.block( i*nr*nr, j, nr*nr, 1) = Eigen::VectorXd::Map(u2.data(), u2.rows()*u2.cols());  
  }  
} 

  return( U );    
}

