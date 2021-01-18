
/*******************************R-callable function***************************************/
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include <iostream>
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/

double mu(double t, arma::mat bb, arma::mat y, double mu0, double theta, arma::mat Q);


double mu(double t, arma::mat bb, arma::mat y, double mu0, double theta, arma::mat Q) {
  arma::mat mu_cmp;
  mu_cmp = (mu0 + bb*y + Q*(y*y))*exp(theta*t);
  return mu_cmp(0,0);
}

RcppExport SEXP complik_discrete(SEXP dat, SEXP n, SEXP m, SEXP mu0, SEXP q, SEXP b, SEXP theta, SEXP k) {
    
    long N = as<long>(n); // Number of rows
    long M = as<long>(m); // Number of columns
    
    arma::mat bb = as<arma::mat>(b);
    arma::mat Q = as<arma::mat>(q);
    double mumu0 = as<double>(mu0); 
    double theta_opt  = as<double>(theta);
    int dim = as<int>(k);
    //Actual data set
    arma::mat dd = as<arma::mat>(dat);  
    double L; // Likelihood
    //End of data loading
    
    arma::mat y1(dim,1);
    arma::mat y2(dim,1);
    
    L = 0;
    for(int i=0; i<N; i++) {
      //Solving differential equations on intervals:
      double t1 = dd(i,1); 
      int jj=0;
      for(int ii=3; ii<M; ii+=2) {
        y1(jj,0) = dd(i,ii);
        y2(jj,0) = dd(i,ii+1);
        jj += 1;
      }
      
      double mu_computed = mu(t1, bb, y1, mumu0, theta_opt, Q);
      
      if(mu_computed < 0) 
        mu_computed = 1e-15;
      if(dd(i,0) == 0) { 
        L += -1*mu_computed;
      } else {
        L += log(1 - exp(-1*mu_computed));
      }
      
    }
    
    
    return(Rcpp::wrap(L));
}

