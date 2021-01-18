// MIT License
// Copyright (c) 2019 Vincent Runge

#ifndef COSTS_H
#define COSTS_H

#include <math.h>

class Costs
{
public:
  Costs();
  double slopeCost(double& u, double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& S2t, double& S2T, double& SPt, double& SPT);
  double vhat(double& v, unsigned int& t, unsigned int& T, double& S1t, double& S1T, double& SPt, double& SPT);
  unsigned int closestStateIndex(double& v, double* states, unsigned int p);
  void fillCoeffsAG(double** coeffs, double* sumY, unsigned int n);
  void linReg(double* coeff, double *gtT, unsigned int t, unsigned int n);
  bool pruningTest(double r, double s, unsigned int tau, unsigned int t, unsigned int n, double S, double Ap, double Am, double Gp, double Gm);
  bool angleTest(unsigned int& t1, unsigned int& t2, unsigned int& t3, double& v1, double& v2, double& v3, double& minAngle);


};
#endif // COSTS_H
