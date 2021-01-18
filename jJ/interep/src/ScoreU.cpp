#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"ScoreU.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace R;

// [[Rcpp::export()]]
Rcpp::List ScoreU(int n, arma::vec& k, arma::vec& y, arma::mat& x, int p, arma::vec& beta, char corre){
  arma::vec aindex=cumsum(k),
    index(aindex.size(),arma::fill::zeros);
  int len=aindex.size()-1;

  for(int i=0; i<len; i++){
    index(i+1)=aindex(i);
  }


  arma::vec U(p,arma::fill::zeros);
  arma::mat dU=arma::mat(p,p,arma::fill::zeros);


    double alfa_hat=0;
    float sum5=0,
      sum6=0,
      sum7=0,
      sum8=0;

    for(int i=0; i<len; i++){
      index(i+1)=aindex(i);
    }
    arma::vec mu=x*beta;
    // vec res=(y-mu)/stddev(y);
    int len1=y.size();
    arma::vec res(len1);

    for(int i=0; i<len1; i++){
      res(i)=(y(i)-mu(i))/stddev(y);
    }

    if(corre == 'i'){
      alfa_hat=0;
    }
    else if(corre == 'e' || corre == 'a'){
      //else if(corre == 'e'){

      for (int i=0; i<n; i++){
        for (int j=0; j<k(i); j++){
          for (int jj=0; jj<k(i); jj++){
            if((j-jj)==1){
              sum7=res(j+index(i))*res(jj+index(i));
              sum5+=sum7;
            }
          }
        }
        sum8=(k(i)-1);
        sum6+=sum8;
      }
      alfa_hat=sum5/sum6;
    }


    int maxclsz=max(k);

    arma::cube Rhat(maxclsz,maxclsz,n);

    for (int i=0; i<n; i++) {
      arma::mat cor1 = arma::mat(k(i),k(i));
      if (corre == 'i') {
        cor1 = arma::mat(k(i),k(i),arma::fill::eye);
      }
      else if (corre == 'e') {
        cor1 = arma::mat(k(i),k(i),arma::fill::ones);
        for(int j=0; j<k(i); j++){
          for(int jj=0; jj<k(i); jj++){
            if(jj != j) {
              cor1(j,jj)=alfa_hat;
            }
          }
        }
      }
      else if (corre == 'a'){
        cor1 = arma::mat(k(i),k(i),arma::fill::ones);
        for(int j=0; j<k(i); j++){
          for(int jj=0; jj<k(i); jj++){
            cor1(j,jj)=pow(alfa_hat,std::abs(j-jj));
          }
        }
      }


      Rhat.slice(i)=cor1;
    }

    for(int i=0; i<n; i++){
      arma::mat A=arma::mat(k(i),k(i),arma::fill::eye);
      arma::vec res(k(i),arma::fill::zeros);
      arma::mat D=arma::mat(k(i),x.n_cols);
      for(int j=0; j<k(i); j++){
        for(int jj=0; jj<p; jj++){
          D(j,jj)=x((j+index(i)),jj);
        }
      }
      arma::vec mu=D*beta;
      for(int t=0; t<k(i); t++){
        res(t)=y(t+index(i))-mu(t);
      }
      arma::mat V=Rhat.slice(i);
      U+=D.t()*inv(V)*res;
      dU+=D.t()*inv(V)*D;
    }

  return Rcpp::List::create(Rcpp::Named("U") = U,
                            Rcpp::Named("dU") = dU);
}

