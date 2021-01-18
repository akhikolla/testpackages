#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List calc_mult_rowsum3(const Eigen::Map<Eigen::ArrayXi> & v, const Eigen::Map<Eigen::ArrayXXd> & B, const Eigen::Map<Eigen::ArrayXXd> & M, const Eigen::Map<Eigen::ArrayXXd> & A, const double ncb2){ 

  //  This function implements:
  //  lapply(1:ncb2, function(i) A * rowsum( B * M  , v))
  //  as before 'v' needs to be sorted!!
  
  const unsigned int n = (unsigned int)ncb2;
  const unsigned int l = v.size();
  const unsigned int mc = M.cols();
  const unsigned int mr = M.rows();
       
  Rcpp::List output(n);

 for (unsigned int u = 0; u < n; u++){
   Eigen::ArrayXXd Res = Eigen::ArrayXXd::Zero(v.maxCoeff() ,mc);
   unsigned int k = 0;
    for (unsigned int i = 0; i < mc; ++i){
      k = 0;
      for (unsigned int j = 0; j < l; ++j){
        Res(k,i) = Res(k,i) + M(j,i) * B(j,u); 
        if(j  < mr - 1 ) {
          if( v[j] != v[j+1] ){ 
            k++;
          } 
        }
      }
    }
    Res = Res.block(0, 0, k + 1, mc);
    Res.array() *= A;
    output[u] = Res;
  }
  return( output );    
}

