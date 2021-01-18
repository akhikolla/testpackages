#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
 
Eigen::MatrixXd calc_mult_rowsum1(const Eigen::Map<Eigen::VectorXi> & v, const Eigen::Map<Eigen::VectorXd> & u, const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::ArrayXXd> & A){         

  //  This function implements:
  //  A *  rowsum( M1 * u , v)
  //  as before 'v' needs to be sorted!!

  const unsigned int l = v.size();
  const unsigned int m = M.cols();  
  const unsigned int mr = M.rows();

  Eigen::MatrixXd Res = Eigen::MatrixXd::Zero(v.maxCoeff() ,m);  
  
  unsigned int k = 0;
  for (unsigned int i = 0; i < m; ++i){
    k = 0;
    for(unsigned int j = 0; j < l; ++j){
      Res(k,i) = Res(k,i) + M(j,i) * u(j); 
      if (j  <  mr - 1 ) {
        if( v[j] != v[j+1] ){ 
          k++;
        } 
      }
    }
  }
  Res = Res.block(0, 0, k + 1, m);
  Res.array() *= A;

  return( Res );    
}

