
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_makeA(const arma::ivec& numSire, const arma::ivec& numDam, const arma::mat& AFounder, const arma::ivec& numFounder, const Rcpp::CharacterVector IndivName){
  int i, j, nSire, nDam, ni, nj;
  int N = numSire.n_elem;
  int NFounder = numFounder.n_elem;
  double tmp;
  Rcpp::NumericMatrix rpedKin(N, N);
  arma::mat pedKin(rpedKin.begin(), rpedKin.nrow(), rpedKin.ncol(), false);
  pedKin.eye();
  
  for(i=0;i<NFounder;i++){
    ni = numFounder.at(i);
    for(j=0;j<=i;j++){
      nj  = numFounder.at(j);
      tmp = AFounder.at(i, j);
      pedKin.at(ni, nj) = tmp;
      pedKin.at(nj, ni) = tmp;
    }
  }
  
  for(i=0;i<N;i++){
    nSire = numSire.at(i);
    nDam  = numDam.at(i);
    if((nSire>0) & (nDam>0)){pedKin.at(i, i) = 1.0 + 0.5*pedKin.at(nSire-1, nDam-1);}
    if((nSire>0) | (nDam>0)){
      for(j=0; j<i; j++){
        tmp = ((nSire>0)?(0.5*pedKin.at(nSire-1, j)):(0.0)) + ((nDam>0)?(0.5*pedKin.at(nDam-1, j)):(0.0));
        pedKin.at(i, j) = tmp;
        pedKin.at(j, i) = tmp;
      }
    }
  }
  
  rpedKin.attr("dimnames") = Rcpp::List::create(IndivName, IndivName);
  return rpedKin;
}
