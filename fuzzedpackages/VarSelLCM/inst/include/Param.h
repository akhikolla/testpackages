/*
Cette classes définie les paramètres pour des données continues

Ces éléments sont:
m_mu : matrice des moyennes
m_sd : matrice des ecarts-types
m_pi : proportions

*/
#ifndef Param_H
#define Param_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Param{
  public:
  Col<double> m_pi;
    
  Param(){};
  ~Param(){};

};
#endif
