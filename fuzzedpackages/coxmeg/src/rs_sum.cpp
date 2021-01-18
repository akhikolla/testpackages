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
Rcpp::List rs_sum(const Eigen::VectorXd & rk_v, const Eigen::VectorXd & d) {
  
  int n = rk_v.size();
  
  Eigen::ArrayXi rs_rs(n);
  Eigen::ArrayXi rs_cs = Eigen::ArrayXi::Zero(n);
  Eigen::ArrayXi rs_cs_p = Eigen::ArrayXi::Zero(n);
  
  for(int i=0; i<n; i++)
  {
    rs_rs(i) = n - rk_v(i);
    rs_cs.tail(rs_rs(i)) += 1;
    if(d(i)>0)
      rs_cs_p.tail(rs_rs(i)) += 1;
  }
  
  return Rcpp::List::create(Rcpp::Named("rs_rs") = rs_rs,
                            Rcpp::Named("rs_cs") = rs_cs,
                            Rcpp::Named("rs_cs_p") = rs_cs_p);
}

