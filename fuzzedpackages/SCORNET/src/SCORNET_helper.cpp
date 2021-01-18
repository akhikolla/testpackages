#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <random>
//#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends("RcppArmadillo")]]


double d_norm(double x, double m, double s){
    double a = (x-m)/s;
    return 1.0/sqrt(2.0*M_PI) * std::exp(-0.5 * a * a);
}


// [[Rcpp::export]]
arma::vec kernelSmoothen(arma::vec hC0, arma::vec C, double b){
    int N = C.n_elem;
    int i,j;
    arma::vec hC0k(N);
    hC0k.zeros();
    
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            hC0k[i] += hC0[j] * d_norm(C[i],C[j],b);
        }
    }
    
    return hC0k;
}

