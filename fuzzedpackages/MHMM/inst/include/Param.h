#ifndef Param_H
#define Param_H

#include "ParamSpecificLaw.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

class Param{
public:
  int m_K, m_M;
  vector< Mat<double> > m_A;
   Col<double> m_delta;
  Mat<double> m_pi;
  ParamSpecificLaw m_lambda;

  Param(){};
  ~Param(){};
  Param(const S4 & obj){
    this->m_K = obj.slot("K");
    this->m_M = obj.slot("M");
    this->m_A = as< vector< Mat<double> >  >(obj.slot("A"));
    this->m_delta = as< vec >(obj.slot("delta"));
    this->m_pi = as< Mat<double> >(obj.slot("pi"));
    this->m_lambda = ParamSpecificLaw(as< List >(obj.slot("lambda")));
  }
};
#endif
