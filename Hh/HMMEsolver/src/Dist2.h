#ifndef _HMMEsolver_DIST2_H
#define _HMMEsolver_DIST2_H


//#include <iostream>
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace arma;

class Dist2{

  int nModel2;
  int n;
  int m;
  int nm;

  double* psi;


public:

  Dist2(int _n, int _m, int _nModel2);
  ~Dist2();


  void Setpsi();
  double* Getpsi();

};


#endif


