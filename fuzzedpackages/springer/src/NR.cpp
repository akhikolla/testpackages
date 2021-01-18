#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"NR.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace R;

// [[Rcpp::export()]]
arma::vec NR(arma::mat& matr1, arma::vec& matr2){

  arma:: vec b1=pinv(matr1)*matr2;




  return b1;
}

