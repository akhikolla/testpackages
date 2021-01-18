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
Rcpp::List ScoreU(int n, arma::vec& k, arma::vec& y, arma::mat& x, int p, arma::vec& beta, char func, char corr){
  arma::vec aindex=cumsum(k),
    index(aindex.size(),arma::fill::zeros);
  int len=aindex.size()-1;

  for(int i=0; i<len; i++){
    index(i+1)=aindex(i);
  }


  arma::vec U(p,arma::fill::zeros);
  arma::mat dU=arma::mat(p,p,arma::fill::zeros);

  if(func=='Q'){
    if(corr=='i'){
      arma::vec phi(p,arma::fill::zeros),
      phi0(p,arma::fill::zeros);
      arma::mat dphi=arma::mat(p,p,arma::fill::zeros),
        C=arma::mat(p,p,arma::fill::zeros);
      for (int i=0; i<n; i++){
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
        phi0=D.t()*res/n;
        phi+=D.t()*res/n;
        dphi-=D.t()*D/n;
        C+=phi0*phi0.t();
      }
      U=2*dphi.t()*pinv(C)*phi;
      dU=2*dphi.t()*pinv(C)*dphi;
    }

    else if(corr=='e'){
      arma::vec phi(2*p,arma::fill::zeros),
      phi0(2*p,arma::fill::zeros),
      phi1(p,arma::fill::zeros),
      phi2(p,arma::fill::zeros);
      arma::mat dphi=arma::mat(2*p,p,arma::fill::zeros),
        dphi0=arma::mat(2*p,p,arma::fill::zeros),
        dphi1=arma::mat(p,p,arma::fill::zeros),
        dphi2=arma::mat(p,p,arma::fill::zeros),
        C=arma::mat(2*p,2*p,arma::fill::zeros);
      for (int i=0; i<n; i++){
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
        arma::mat M1=arma::mat(k(i),k(i),arma::fill::ones);
        for (int u=0;u<k(i);u++){
          for (int v=0;v<k(i);v++){
            if(u==v) M1(u,v)=0;
          }
        }
        phi1=D.t()*res/n;
        phi2=D.t()*M1*res/n;
        dphi1=D.t()*D/n;
        dphi2=D.t()*M1*D/n;
        for (int u=0;u<p;u++){
          phi0(u)=phi1(u);
          for (int v=0;v<p;v++){
            dphi0(u,v)=dphi1(u,v);
          }
        }
        for (int u=p;u<2*p;u++){
          phi0(u)=phi2(u-p);
          for (int v=0;v<p;v++){
            dphi0(u,v)=dphi2(u-p,v);
          }
        }
        phi+=phi0;
        dphi-=dphi0;
        C+=phi0*phi0.t();
      }
      //std::cout << C(0,0) << endl;
      U=2*dphi.t()*pinv(C)*phi;
      dU=2*dphi.t()*pinv(C)*dphi;
    }

    else if(corr=='a'){
      arma::vec phi(3*p,arma::fill::zeros),
      phi0(3*p,arma::fill::zeros),
      phi1(p,arma::fill::zeros),
      phi2(p,arma::fill::zeros),
      phi3(p,arma::fill::zeros);
      arma::mat dphi=arma::mat(3*p,p,arma::fill::zeros),
        dphi0=arma::mat(3*p,p,arma::fill::zeros),
        dphi1=arma::mat(p,p,arma::fill::zeros),
        dphi2=arma::mat(p,p,arma::fill::zeros),
        dphi3=arma::mat(p,p,arma::fill::zeros),
        C=arma::mat(3*p,3*p,arma::fill::zeros);
      for (int i=0; i<n; i++){
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
        arma::mat M1=arma::mat(k(i),k(i),arma::fill::zeros),
          M2=arma::mat(k(i),k(i),arma::fill::zeros);
        for (int u=0;u<k(i);u++){
          for (int v=0;v<k(i);v++){
            if((u-v)==1 || (v-u)==1) M1(u,v)=1;
          }
        }
        M2(0,0)=1;
        M2(k(i)-1,k(i)-1)=1;
        //std::cout << M1 << endl;

        phi1=D.t()*res/n;
        phi2=D.t()*M1*res/n;
        phi3=D.t()*M2*res/n;
        dphi1=D.t()*D/n;
        dphi2=D.t()*M1*D/n;
        dphi3=D.t()*M2*D/n;
        for (int u=0;u<p;u++){
          phi0(u)=phi1(u);
          for (int v=0;v<p;v++){
            dphi0(u,v)=dphi1(u,v);
          }
        }
        for (int u=p;u<(2*p);u++){
          phi0(u)=phi2(u-p);
          for (int v=0;v<p;v++){
            dphi0(u,v)=dphi2(u-p,v);
          }
        }
        for (int u=(2*p);u<(3*p);u++){
          phi0(u)=phi3(u-(2*p));
          for (int v=0;v<p;v++){
            dphi0(u,v)=dphi3(u-(2*p),v);
          }
        }
        phi+=phi0;
        dphi-=dphi0;
        C+=phi0*phi0.t();
      }

      U=2*dphi.t()*pinv(C)*phi;
      dU=2*dphi.t()*pinv(C)*dphi;
    }
  }

  else if(func=='G'){
    double alfa_hat=0;
    float sum5=0,
      sum6=0,
      sum7=0,
      sum8=0;
    int len=aindex.size()-1;
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

    if(corr == 'i'){
      alfa_hat=0;
    }
    else if(corr == 'e' || corr == 'a'){
      //else if(corr == 'e'){

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
      if (corr == 'i') {
        cor1 = arma::mat(k(i),k(i),arma::fill::eye);
      }
      else if (corr == 'e') {
        cor1 = arma::mat(k(i),k(i),arma::fill::ones);
        for(int j=0; j<k(i); j++){
          for(int jj=0; jj<k(i); jj++){
            if(jj != j) {
              cor1(j,jj)=alfa_hat;
            }
          }
        }
      }
      else if (corr == 'a'){
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
  }
  return Rcpp::List::create(Rcpp::Named("U") = U,
                            Rcpp::Named("dU") = dU);
}

