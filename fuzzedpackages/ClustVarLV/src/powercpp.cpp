#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]

Rcpp::List powerEigen(Eigen::MatrixXd & X) {
    double temp=0.;
    double change=1. ;
    double lambda=0.;
    int niter=0 ;
    Eigen::VectorXd v= Eigen::VectorXd::Constant(X.cols(),1);
    while( change > 0.000001 && niter < 1000 ){
        v = X*v;
        v /= v.norm() ;
        lambda = v.dot( X*v);
        if( niter == 0)
             temp = lambda;
        else if( temp != 0.){
             change = fabs (1. - lambda/temp );
             temp = lambda;
        }
        niter ++;
   }

   return Rcpp::List::create(Rcpp::Named("values")  = lambda,
                             Rcpp::Named("vectors") = v       );
 }
