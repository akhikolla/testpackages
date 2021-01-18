
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd calc_rowsum(const Eigen::Map<Eigen::VectorXi> & v , const Eigen::Map<Eigen::MatrixXd> & M){ 
        
  //  This function implements:
  //  the computation of column sums across rows of a matrix M for each level/integer number of the variable v. 
  //  this function requires 'v' to be sorted and dense. DO NOT USE GENERICALLY.

  const unsigned int l = v.size();
  const unsigned int m = M.cols();  

  Eigen::MatrixXd Res = Eigen::MatrixXd::Zero(v.maxCoeff() ,m);  
  //Rcpp::Rcout << "Got here!\n"; 
  // for each column
  unsigned int k = 0;
  for (unsigned int i = 0; i < m; ++i){
    k = 0;
    for(unsigned int j = 0; j < l; ++j){
      // Rcpp::Rcout << " i, j, k: " << i<< ", " << j << ", " << k << " \n";
      Res(k,i) = Res(k,i) + M(j,i); 
      if (j < M.rows()-1){
        if( v[j] != v[j+1] ){ 
          k++;
        }
      }
    }
  }
  return( Res.block(0, 0, k + 1, m) );    
}

