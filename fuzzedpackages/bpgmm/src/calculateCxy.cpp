// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List Calculate_Cxy(int m,
                         int n,
                         Rcpp::S4 hparam,
                         Rcpp::S4 thetaYList,
                         arma::vec ZOneDim,
                         arma::vec qVec,
                         arma::mat X){

  double alpha1 = hparam.slot("alpha1");
  double alpha2 = hparam.slot("alpha2");

  List Y      = thetaYList.slot("Y");
  List lambda = thetaYList.slot("lambda");
  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");

  arma::mat Zmat = get_Z_mat(ZOneDim, m, n);
  List A(m);
  arma::vec nVec(m);

  for(int k=0; k<m; ++k) {
    nVec(k) = sum(Zmat.row(k));
    arma::vec alpha1_vec(1);

    A(k) =  diagmat(join_cols(alpha1_vec.fill(alpha1), arma::ones(qVec(k)) * alpha2));
  };

  List Cxxk;
  List Cxyk;
  List Cyyk;
  List Cytytk;
  List Cxtytk;
  List CxL1k;
  List Cxmyk;

  for (int k=0; k<m; ++k) {
    arma::mat Cxxkk;
    arma::mat Cxykk;
    arma::mat Cyykk;
    arma::mat Cytytkk;
    arma::mat Cxtytkk;
    arma::mat CxL1kk;
    arma::mat Cxmykk;
    arma::mat y_k = Y(k);
    arma::vec m_k = M(k);
    arma::mat lambda_k = lambda(k);
    arma::vec one_vec(1);

    for (int i=0; i<n; ++i) {
      if(i == 0){
        Cxxkk =  Zmat(k,i) * (X.col(i) * trans(X.col(i)));
        Cxykk =  Zmat(k,i) * (X.col(i) * trans(y_k.col(i)));
        Cyykk =  Zmat(k,i) * (y_k.col(i) * trans(y_k.col(i)));
        Cxtytkk =  Zmat(k,i) * (X.col(i) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cytytkk =  Zmat(k,i) * (join_cols(one_vec.fill(1), y_k.col(i)) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cxmykk =  Zmat(k,i) * ((X.col(i) - m_k) * trans(y_k.col(i)));
        CxL1kk = Zmat(k,i) * (X.col(i) -  (lambda_k * y_k.col(i)));

      }else{
        Cxxkk = Cxxkk + Zmat(k,i) * (X.col(i) * trans(X.col(i)));
        Cxykk = Cxykk + Zmat(k,i) * (X.col(i) * trans(y_k.col(i)));
        Cyykk = Cyykk +  Zmat(k,i) * (y_k.col(i) * trans(y_k.col(i)));
        Cxtytkk = Cxtytkk + Zmat(k,i) * (X.col(i) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cytytkk = Cytytkk + Zmat(k,i) * (join_cols(one_vec.fill(1), y_k.col(i)) * trans(join_cols(one_vec.fill(1), y_k.col(i))));
        Cxmykk = Cxmykk + Zmat(k,i) * ((X.col(i) - m_k) * trans(y_k.col(i)));
        CxL1kk = CxL1kk + Zmat(k,i) * (X.col(i) -  (lambda_k * y_k.col(i)));
      }
    }
    Cxxk.push_back(Cxxkk);
    Cxyk.push_back(Cxykk);
    Cyyk.push_back(Cyykk);
    Cxtytk.push_back(Cxtytkk);
    Cytytk.push_back(Cytytkk);
    Cxmyk.push_back(Cxmykk);
    CxL1k.push_back(CxL1kk);
  }

  arma::mat sumCxmyk;
  arma::mat sumCyyk;

  for (int k=0; k<m; ++k) {
      arma::mat cxmyk_k = Cxmyk[k];
      arma::mat cyyk_k = Cyyk[k];
    if (k == 0) {
      sumCxmyk = cxmyk_k;
      sumCyyk  = cyyk_k + diagmat( arma::ones(qVec(k))) * alpha2;
    } else {
      sumCxmyk = sumCxmyk + cxmyk_k;
      sumCyyk  = sumCyyk + cyyk_k + diagmat(arma::ones(qVec(k))) * alpha2;
    }
  }

  List res = Rcpp::List::create(Named("A")        = A,
                                Named("nVec")     = nVec,
                                Named("Cxxk")     = Cxxk,
                                Named("Cxyk")     = Cxyk,
                                Named("Cyyk")     = Cyyk,
                                Named("Cytytk")   = Cytytk,
                                Named("Cxtytk")   = Cxtytk,
                                Named("CxL1k")    = CxL1k,
                                Named("Cxmyk")    = Cxmyk,
                                Named("sumCxmyk") = sumCxmyk,
                                Named("sumCyyk")  = sumCyyk);

  return(res);
}
