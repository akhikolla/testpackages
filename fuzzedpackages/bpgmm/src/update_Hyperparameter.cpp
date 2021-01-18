// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
S4 update_Hyperparameter(
    int m,
    int p,
    int q,
    Rcpp::S4 hparam,
    Rcpp::S4 thetaYList,
    arma::vec dVec,
    arma::vec sVec
  ){


  // extract value from s4 obj
  double delta   = hparam.slot("delta");
  double ggamma  = hparam.slot("ggamma");
  List M         = thetaYList.slot("M");
  List lambda    = thetaYList.slot("lambda");
  List psy       = thetaYList.slot("psy");

  // update alpha1
  double alpha1Rate = 0;
  arma::vec alpha1RateTemp;
  for(int k = 0; k < m; k++){
    arma::vec Mk = M(k);
    arma::mat psyk = psy(k);
    alpha1RateTemp = 0.5 * trans(Mk) * arma::inv(psyk) * Mk;
    alpha1Rate += alpha1RateTemp(0);
  }
  alpha1Rate += sVec(0);

  double alpha1Shape = m*p/2.0 + dVec(0);
  double alpha1Scale = 1.0/alpha1Rate;

  // Rcpp::rgamma(n, shape  , scale), scale = 1/rate;
  arma::vec alpha1Vec =  Rcpp::rgamma(1, alpha1Shape, alpha1Scale);

  double alpha2Rate = 0;
  arma::mat alpha2RateMat = arma::zeros(q,q);
  arma::mat alpha2RateTemp;
  for(int k = 0; k < m; k++){

    arma::mat lambdak = lambda(k);
    arma::mat psyk = psy(k);
    alpha2RateTemp = 0.5 * trans(lambdak) * arma::inv(psyk)  * lambdak;
    alpha2RateMat += alpha2RateTemp;
  }

  alpha2Rate += sVec(1) + sum(alpha2RateMat.diag());

  double alpha2Shape = q*p/2.0 + dVec(1);
  double alpha2Scale = 1.0/alpha2Rate;

  //Rcpp::rgamma(n, shape  , scale), scale = 1/rate;
  arma::vec alpha2Vec =  Rcpp::rgamma(1, alpha2Shape, alpha2Scale);


  double bbetaRate = 0;
  arma::mat bbetaRateMat = arma::zeros(p,p);
  arma::mat bbetaRateTemp;
  for(int k = 0; k < m; k++){
    arma::mat psyk = psy(k);
    bbetaRateTemp = arma::inv(psyk);
    bbetaRate  += sum(bbetaRateTemp.diag());
  }
  bbetaRate += dVec(2);

  double bbetaShape = m*p*delta + dVec(2);
  double bbetaScale = 1.0/bbetaRate;

  arma::vec bbetaVec =  Rcpp::rgamma(1, bbetaShape, bbetaScale);
  // Creating an object of Hparam class
  S4 newhparam("Hparam");

  // Setting values to the slots
  newhparam.slot("alpha1")  = alpha1Vec(0);
  newhparam.slot("alpha2")  = alpha2Vec(0);
  newhparam.slot("delta")   = delta;
  newhparam.slot("ggamma")  = ggamma;
  newhparam.slot("bbeta")   = bbetaVec(0);

  return(newhparam);
}
