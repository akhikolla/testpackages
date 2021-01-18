#ifndef STCparam_H
#define STCparam_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class STCparam{
  public:
    Mat<double> m_sigma;
    vector< Mat<double> > m_lambda, m_beta;
    Col<double> m_proportions;
  
  STCparam(){};
  STCparam(const S4 & obj){
    this->m_proportions = as<vec>(obj.slot("proportions"));
    this->m_sigma  = as<mat>(obj.slot("sigma")); 
    this->m_lambda.resize(m_proportions.n_rows);
    this->m_beta.resize(m_proportions.n_rows);
    List tmplambda = List(obj.slot("lambda"));
    List tmpbeta = List(obj.slot("beta"));
    for (int g=0; g<m_proportions.n_rows; g++){
      this->m_lambda[g] = as<mat>(tmplambda[g]);
      this->m_beta[g] = as<mat>(tmpbeta[g]);
    }
  };
  ~STCparam(){};  
};
#endif
