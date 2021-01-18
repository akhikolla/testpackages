// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::IntegerVector update_PostZ(
                     arma::mat X,
                     int m,
                     int n,
                     Rcpp::S4 thetaYList){

  List lambda   = thetaYList.slot("lambda");
  List Y        = thetaYList.slot("Y");
  List M        = thetaYList.slot("M");
  List psy      = thetaYList.slot("psy");
  arma::vec tao = thetaYList.slot("tao");




  arma::mat pMat(m, n);
  arma::mat dMat(m, n);

 //   std::cout  << k << std::endl;


//std::cout  << "1" << std::endl;
  for(int k = 0; k < m; k++ ){
    for(int i = 0; i < n; i++){

      arma::vec Mk = M(k);
      arma::mat lambdak = lambda(k);
      arma::mat psyk = psy(k);
      arma::mat var = psyk + lambdak * trans(lambdak);
       // std::cout  << k << std::endl;
      arma::vec temp = dmvnrm_arma(trans(X.col(i)), trans(Mk), var, true);
      dMat(k, i) = temp.at(0);
    // try {
    //   arma::vec temp = dmvnrm_arma(trans(X.col(i)), trans(Mk), var, true);
    //   dMat(k, i) = temp.at(0);
    // } catch(const std::exception & e){
    //     dMat(k, i) = -100000;
    // }


    //  dMat(k, i) = temp.at(0);
    }
  }
    // std::cout  << dMat(0, 0) << std::endl;
//     std::cout  << dMat(1, 0) << std::endl;
//     std::cout  << dMat(2, 0) << std::endl;
//     std::cout  << dMat(3, 0) << std::endl;
//   std::cout  << "2" << std::endl;
  for(int k = 0; k < m; k++ ){
    double taok = tao(k);
    arma::vec taokvec = rep(taok, n);
    dMat.row(k) += trans(taokvec);
  }


    for(int i = 0; i < n; i++){
      for(int k = 0; k < m; k++ ){
        pMat(k, i) = calculate_Ratio(dMat(k, i),dMat.col(i));
      }
    }

  Rcpp::IntegerVector ZOneDim(n);
  Rcpp::IntegerVector ZOneDimTemp;
  arma::vec tempProb;
  Rcpp::IntegerVector sampleVec = Rcpp::seq(1, m);

    for(int i = 0; i < n; i++){
      tempProb =  pMat.col(i);
      ZOneDimTemp = RcppArmadillo::sample(sampleVec, 1, FALSE, tempProb);
      ZOneDim(i) = ZOneDimTemp(0);
    }

  return(ZOneDim);
}
