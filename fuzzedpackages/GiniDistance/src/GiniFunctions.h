#ifndef _INCL_GiniFunctions_
#define _INCL_GiniFunctions_
#include <RcppArmadillo.h>
using namespace Rcpp;

// sum of vector
double VectorSum(NumericVector x);

// variance of vector
double VectorVar(NumericVector x);

// order a vector
IntegerVector orderc(NumericVector x);

//subset a matrix
NumericMatrix ss(NumericMatrix X_, IntegerVector ind_);


// Rcpp distance matrix
NumericMatrix rcpp_Eu_distance(NumericMatrix mat);

// Rcpp Kernel Distance matrix
NumericMatrix rcpp_Kernel_Distance(NumericVector vec, double sigma);

// Rcpp vector GiniDistance 
double rcpp_covg(NumericVector x, NumericVector y);

// Rcpp vector Gini Covariance to power alpha 
double rcpp_gCov_alpha(NumericVector x, NumericVector y, double alpha);

// Rcpp vector Gini Correlation to power alpha 
double rcpp_gCor_alpha(NumericVector x, NumericVector y, double alpha);

// Rcpp vector GiniDistance power to alpha
double rcpp_covg_alpha(NumericVector x, NumericVector y, double alpha);

// Rcpp vector Kernel GiniDistance power to alpha
double rcpp_Kcovg_alpha(NumericVector x, NumericVector y, double sigma, double alpha);

// Rcpp matrix GiniDistance  
double Rcpp_Covg(NumericMatrix x, NumericVector y);

// Rcpp matrix GiniDistance raise to power alpha 
double Rcpp_Covg_Alpha(NumericMatrix x, NumericVector y, double alpha);

// Rcpp matrix Kernel GiniDistance raise to power alpha 
double Rcpp_KCovg_Alpha(NumericMatrix x, NumericVector y, double sigma, double alpha);

// Rcpp matrix Kernel GiniDistance multi classes raise to power alpha 
double Rcpp_KgCov_Alpha(NumericMatrix x,double sigma, double alpha);

// Rcpp matrix Kernel Gini Correlation multi classes raise to power alpha 
double Rcpp_KgCor_Alpha(NumericMatrix x,double sigma, double alpha);

// Rcpp matrix GiniCovariance multi-classes
double Rcpp_gCov(NumericMatrix x); 

// Rcpp matrix GiniCovariance multi-classes power alpha
double Rcpp_gCov_Alpha(NumericMatrix x, double alpha); 

// Rcpp matrix Ginicorrelation multi-classes power alpha
double Rcpp_gCor_Alpha(NumericMatrix x, double alpha); 


// Rcpp matrix GiniCorrelation multi-classes
double Rcpp_gCor(NumericMatrix x); 

//RcppKernelGiniDistance for vector
double rcpp_Kcovg(NumericVector x, NumericVector y, double sigma);

//RcppKernelGiniDistance for matrix
double Rcpp_KCovg(NumericMatrix x, NumericVector y, double sigma);

//RcppKernelGiniCovariance
double Rcpp_KgCov(NumericMatrix x, double sigma);

//RcppGiniCovariance for two vectors
double rcpp_gCov(NumericVector x, NumericVector y);

//RcppGiniCorrelation for two vectors
double rcpp_gCor(NumericVector x, NumericVector y);

//RcppKernelGiniCorrelation
double Rcpp_KgCor(NumericMatrix x, double sigma);

//jackknife variance  
double Rcpp_HatV_gCov(NumericMatrix x);

//jackknife variance to power alpha  
double Rcpp_HatV_gCov_Alpha(NumericMatrix x, double alpha);

//jackknife variance of correlation to power alpha  
double Rcpp_HatV_gCor_Alpha(NumericMatrix x, double alpha);

//jackknife variance  
double Rcpp_HatV_gCor(NumericMatrix x);

//jackknife variance of Kernel Gini covariance  
double Rcpp_HatV_KgCov(NumericMatrix x, double sigma);

//jackknife variance of Kernel Gini covariance to power alpha 
double Rcpp_HatV_KgCov_Alpha(NumericMatrix x, double sigma, double alpha);

//jackknife variance of Kernel Gini correlation  
double Rcpp_HatV_KgCor(NumericMatrix x, double sigma);


#endif

