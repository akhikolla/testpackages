
#include <RcppEigen.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd fast_lapply_length(Rcpp::List const input1, Rcpp::List const input2, Rcpp::NumericVector const Ind){ 
          
  //  This function implements:
  //  do.call(rbind, lapply(Ind, function(i) input1[[i]] %*% input2[[i]] ))  

  const unsigned int l = Ind.size();
  unsigned int s = 0; 

  for (unsigned int i = 0; i != l; ++i){  
    Rcpp::NumericMatrix xx = input1[ Ind(i) ];
    s += xx.nrow();
  }

  Rcpp::NumericMatrix xx = input2[1];

  const unsigned int nc = xx.ncol();
  Eigen::MatrixXd U( s , nc );
 
  unsigned int j = 0;
  for (unsigned int i = 0; i != l; ++i){  
    const Eigen::MatrixXd M1 = (input1[ Ind(i) ]);
    const Eigen::MatrixXd M2 = (input2[ Ind(i) ]);
    U.block( j, 0, M1.rows() , nc ) =   M1 * M2;
    j += M1.rows();
  }
  return( U );   
}

