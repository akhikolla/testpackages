#ifndef STCtune_H
#define STCtune_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class STCtune{
  public:
    int m_nbinitSmall, m_nbinitKept, m_nbiterSmall, m_nbiterKept;
    double m_tol;
    
  STCtune(){};
  STCtune(const S4 & obj){
    this->m_tol = obj.slot("tol");
    this->m_nbinitSmall  = obj.slot("nbinitSmall"); 
    this->m_nbinitKept= obj.slot("nbinitKept");
    this->m_nbiterSmall = obj.slot("nbiterSmall");
    this->m_nbiterKept = obj.slot("nbiterKept");
  };
  ~STCtune(){};  
};
#endif
