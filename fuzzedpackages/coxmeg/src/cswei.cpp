// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Eigen::VectorXd cswei(const Eigen::Map<Eigen::VectorXd> w_v, const Eigen::Map<Eigen::VectorXd> rs_rs, 
                      const Eigen::MatrixXi & ind, const Eigen::VectorXd & rev) {
  
  int n = w_v.size();
  
  Eigen::VectorXd temp(n);
  for(int j=0; j<n; j++)
    temp(j) = w_v(ind(j,0));
  if(rev(0)>0)
  {
    temp = temp.reverse().eval();
  }
  for(int j=1; j<n; ++j)
    temp(j) += temp(j-1);
  Eigen::VectorXd w_v_x(n);
  for(int j=0; j<n; j++)
    w_v_x(j) = temp(rs_rs(j));
  for(int j=0; j<n; j++)
    temp(j) = w_v_x(ind(j,1));
  
  return temp;
  
}

