#ifndef STCmodel_H
#define STCmodel_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class STCmodel{
  public:
    int m_G, m_K, m_Q, m_spatial, m_nbparam;
  
  STCmodel(){};
  STCmodel(const S4 & obj){
    this->m_G = obj.slot("G");
    this->m_K  = obj.slot("K"); 
    this->m_Q= obj.slot("Q");
    this->m_spatial= obj.slot("spatial");
    this->m_nbparam = obj.slot("nbparam");
  };
  ~STCmodel(){};  
};
#endif
