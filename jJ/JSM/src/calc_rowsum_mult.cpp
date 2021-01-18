#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd calc_rowsum_mult(const Eigen::Map<Eigen::VectorXi> & v, const Eigen::Map<Eigen::VectorXd> & u, const Eigen::Map<Eigen::MatrixXd> & M){ 
     

  //  This function implements:
  //  rowsum( M1 * u , v)
  //  as before 'v' needs to be sorted!!
    
  const unsigned int l = v.size();
  const unsigned int m = M.cols();  

  if ( l != u.size() ) {  
    Rcpp::stop("The multiplier vector u and the grouping vector v need to have the same size.\n");
  }
  Eigen::MatrixXd Res = Eigen::MatrixXd::Zero(v.maxCoeff() ,m);  
  
  unsigned int k = 0;
  for (unsigned int i = 0; i < m; ++i){
    k = 0;
    for(unsigned int j = 0; j < l; ++j){
 //     Rcpp::Rcout << "i, j, k: " << i << ", " << j << ", " << k <<"\n";
      Res(k,i) = Res(k,i) + M(j,i) * u(j); 
      if (j  < M.rows() - 1 ) {
        if( v[j] != v[j+1] ){ 
          k++;
        } 
      }
    }
  }

  return( Res.block(0, 0, k + 1, m) );    
}

