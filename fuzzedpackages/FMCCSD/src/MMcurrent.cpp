#include <RcppArmadillo.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 





// [[Rcpp::export]]
double likelihoodi(const double&b,const arma::vec&parameters,
                              const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                              const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  double result=0;double S=0;arma::mat midresult(1,1);
  arma::vec covariate(betadim+gammadim);
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    if(r!=0){
      S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                    *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
    }
      else{
      S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
      }

    if(S>0.99999999999){
      S=0.99999999999;
    }

    if((S)<std::pow(10,-30)){
      
      S=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
  }
  
    result=std::exp(result);
  // result=result*std::exp(b*b);
  return result;
}
// [[Rcpp::export]]
arma::mat weightfunction(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                         const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                         const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim){
  int order=rules.n_rows;
  arma::vec normvec=zeros(order);
  for(int i=0;i<order;i++){
    normvec(i)=R::dnorm(rules(i,0),0,1,0);
  }
  arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
  arma::mat result=zeros(order,n);
  arma::vec functionvalue(order);
  for(int i=0;i<n;i++){
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodi(rules(k,0),parameters,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
      
    }
    result.col(i)=weightvec%functionvalue/sum(weightvec%functionvalue);
  }
  return result;
}
  
// [[Rcpp::export]]
double MMfunctioni(const double&b,const arma::vec&parameters,const arma::vec&lastpar,
                   const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                   const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  arma::vec firstpart=zeros(ni);
  arma::vec secondpart=zeros(ni);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,parameters.n_elem-1);
  arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,parameters.n_elem-1);
  arma::vec estbeta=parameters.subvec(0,betadim-1);
  arma::vec lastbeta=lastpar.subvec(0,betadim-1);
  arma::vec estgamma=parameters.subvec(betadim,betadim+gammadim-1);
  arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
  double esttheta=parameters(betadim+gammadim);
  double lasttheta=lastpar(betadim+gammadim);
  arma::vec estpara=zeros(ni);
  arma::vec lastpara=zeros(ni);
  double result=0;
  if(r!=0){
  estpara=1+r*(trans(blC)*exp(estpsi))%exp(X*estbeta+sum(trans(estgamma)*Z)+exp(esttheta)*b);
  lastpara=1+r*(trans(blC)*exp(lastpsi))%exp(X*lastbeta+sum(trans(lastgamma)*Z)+exp(lasttheta)*b);
  firstpart=Delta%(1-(1-pow(lastpara,-1/r))%pow(1-pow(estpara,-1/r),-1));
  secondpart=(1/r)*(1-Delta)%(1-estpara%pow(lastpara,-1));
  result=sum(firstpart+secondpart);
  
  
  }
  else{
    estpara=(trans(blC)*exp(estpsi))%exp(X*estbeta+sum(trans(estgamma)*Z)+exp(esttheta)*b);
    lastpara=(trans(blC)*exp(lastpsi))%exp(X*lastbeta+sum(trans(lastgamma)*Z)+exp(lasttheta)*b);
    firstpart=Delta%(1-(1-exp(-lastpara))%pow((1-exp(-estpara)),-1));
    secondpart=(1-Delta)%(lastpara-estpara);
    result=sum(firstpart+secondpart);
  }
  return result;
}
  
// [[Rcpp::export]]
double targetfunc(const arma::vec&parameters,const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&C,const arma::field<arma::vec>&Delta,
                  const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                  const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                  const arma::mat&R,const double&lambda){
  int order=rules.n_rows;
  arma::mat weightmat=weightfunction(lastpar,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim);
  arma::mat valuemat=zeros(order,n);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,parameters.n_elem-1);
  for(int i=0;i<n;i++){
    for(int k=0;k<order;k++){
      valuemat(k,i)=MMfunctioni(rules(k,0),parameters,lastpar,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
  }
  
  double result=accu(weightmat%valuemat);
  
  
    arma::mat penalty(estpsi.n_elem,estpsi.n_elem);
    penalty=exp(trans(estpsi))*R*exp(estpsi);
    result=result-lambda*penalty(0,0);
  
  return -result;
}


// [[Rcpp::export]]
double MMfunctioninofrailty(const double&b,const arma::vec&parameters,const arma::vec&lastpar,
                   const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                   const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  arma::vec firstpart=zeros(ni);
  arma::vec secondpart=zeros(ni);
  arma::vec estpsi=parameters.subvec(betadim+gammadim,parameters.n_elem-1);
  arma::vec lastpsi=lastpar.subvec(betadim+gammadim,parameters.n_elem-1);
  arma::vec estbeta=parameters.subvec(0,betadim-1);
  arma::vec lastbeta=lastpar.subvec(0,betadim-1);
  arma::vec estgamma=parameters.subvec(betadim,betadim+gammadim-1);
  arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
  arma::vec estpara=zeros(ni);
  arma::vec lastpara=zeros(ni);
  double result=0;
  if(r!=0){
    estpara=1+r*(trans(blC)*exp(estpsi))%exp(X*estbeta+sum(trans(estgamma)*Z));
    lastpara=1+r*(trans(blC)*exp(lastpsi))%exp(X*lastbeta+sum(trans(lastgamma)*Z));
    firstpart=Delta%(1-(1-pow(lastpara,-1/r))%pow(1-pow(estpara,-1/r),-1));
    secondpart=(1/r)*(1-Delta)%(1-estpara%pow(lastpara,-1));
    result=sum(firstpart+secondpart);
    
    
  }
  else{
    estpara=(trans(blC)*exp(estpsi))%exp(X*estbeta+sum(trans(estgamma)*Z));
    lastpara=(trans(blC)*exp(lastpsi))%exp(X*lastbeta+sum(trans(lastgamma)*Z));
    firstpart=Delta%(1-(1-exp(-lastpara))%pow((1-exp(-estpara)),-1));
    secondpart=(1-Delta)%(lastpara-estpara);
    result=sum(firstpart+secondpart);
  }
  return result;
}

// [[Rcpp::export]]
double targetfuncfrailty(const arma::vec&parameters,const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&C,const arma::field<arma::vec>&Delta,
                  const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                  const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                  const arma::mat&R,const double&lambda){
  int order=rules.n_rows;
  arma::mat weightmat=weightfunction(lastpar,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim);
  arma::mat valuemat=zeros(order,n);
  arma::vec estpsi=parameters.subvec(betadim+gammadim,parameters.n_elem-1);
  for(int i=0;i<n;i++){
    for(int k=0;k<order;k++){
      valuemat(k,i)=MMfunctioninofrailty(rules(k,0),parameters,lastpar,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
  }
  
  double result=accu(weightmat%valuemat);
  
    arma::mat penalty(estpsi.n_elem,estpsi.n_elem);
    penalty=exp(trans(estpsi))*R*exp(estpsi);
    result=result-lambda*penalty(0,0);
  
  return -result;
}

// [[Rcpp::export]]
double likelihoodfunc1current(const double&b,const arma::vec&parameters,const arma::vec&C,
                              const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                              const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  double result=0;double S=0;
  arma::mat midresult(1,1);
  arma::vec covariate(betadim+gammadim);
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
    
    
    if(r==0){
      S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
    }
    else{
      S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                    *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
    }
    
    
    if(S>0.99999999999){
      S=0.99999999999;
    }
    if((S)<std::pow(10,-30)){
      
      S=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
  }
  
  
  result=result+R::dnorm(b,0,1,true);
  result=std::exp(result);
  result=result*std::exp(b*b);
  return result;
}

// This function is used to caculate the likelihood value using Hermit quadrature.

// [[Rcpp::export]]
double testquadrature1current(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&C,const arma::field<arma::vec>&Delta,
                              const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                              const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                              const arma::mat&R,const double&lambda){
  int totaldim=parameters.n_elem;
  int order=rules.n_rows;double result=0;
  arma::vec weightvec;weightvec=rules.col(1);double term1=0;
  arma::vec functionvalue(order);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,totaldim-1);
  for(int i=0;i<n;i++){
    
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodfunc1current(rules(k,0),parameters,C(i),
                    Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
    
    term1=sum(functionvalue%weightvec);
    if(term1<std::pow(10,-30)){
      term1=std::pow(10,-30);}
    result=result+std::log(term1);
    
  }
  
    arma::mat penalty(estpsi.n_elem,estpsi.n_elem);
    penalty=exp(trans(estpsi))*R*exp(estpsi);
    result=result-lambda*penalty(0,0);
  
  return -result;
}

// [[Rcpp::export]]
double penaltyterm(const arma::vec&psi,const double&lambda,const arma::mat&R){
  arma::mat penalty=exp(trans(psi))*R*exp(psi);
  double result=lambda*penalty(0,0);
  return(result);
}







// [[Rcpp::export]]
double likelihoodfunc1currentnofrailty(const double&b,const arma::vec&parameters,const arma::vec&C,
                                       const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                                       const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
  int totaldim=parameters.n_elem;
  double result=0;double S=0;arma::mat midresult;midresult.zeros(1,1);
  arma::vec covariate(betadim+gammadim);
  for(int j=0;j<ni;j++){
    covariate.subvec(0,betadim-1)=trans(X.row(j));
    covariate.subvec(betadim,betadim+gammadim-1)=Z;
    
    // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
    
    
    if(r==0){
      S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                   *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)));
    }
    else{
      S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                    *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate))),-1/r);
    }
    
    
    if(S>0.99999999999){
      S=0.99999999999;
    }

    if((S)<std::pow(10,-30)){
      
      S=std::pow(10,-30);
    }
    result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
  }
  
  
  result=result+R::dnorm(b,0,1,true);
  result=std::exp(result);
  result=result*std::exp(b*b);
  return result;
}

// This function is used to caculate the likelihood value using Hermit quadrature.

// [[Rcpp::export]]
double testquadrature1currentnofrailty(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&C,const arma::field<arma::vec>&Delta,
                                       const arma::field<arma::vec>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                                       const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                                       const arma::mat&R,const double&lambda){
  int totaldim=parameters.n_elem;
  int order=rules.n_rows;double result=0;
  arma::vec weightvec;weightvec=rules.col(1);double term1=0;
  arma::vec functionvalue(order);
  arma::vec estpsi=parameters.subvec(betadim+gammadim+1,totaldim-1);
  
  for(int i=0;i<n;i++){
    
    for(int k=0;k<order;k++){
      functionvalue(k)=likelihoodfunc1currentnofrailty(rules(k,0),parameters,C(i),
                    Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
    }
    
    term1=sum(functionvalue%weightvec);
    if(term1<std::pow(10,-30)){
      term1=std::pow(10,-30);}
    result=result+std::log(term1);
    
  }
  arma::mat penalty(estpsi.n_elem,estpsi.n_elem);
  penalty=exp(trans(estpsi))*R*exp(estpsi);
  result=result-lambda*penalty(0,0);
  
  return -result;
}








