#ifndef _HMMEsolver_DIST_H
#define _HMMEsolver_DIST_H


//#include <iostream>
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace arma;

class Dist{

  int nModel;
  int n;
  int m;
  int nm;

  double* W;
  double* eta;
  double* deta_dmu;
  double* s;


public:

  Dist(int _n, int _m, int _nModel);
  ~Dist();

  void SetW(double* _mu);
  void Seteta(double* _mu);
  void Setdeta_dmu(double* _mu);
  void Sets(double* _mu, double* _y);


  double* GetW();
  double* Geteta();
  double* Getdeta_dmu();
  double* Gets();

};


#endif
