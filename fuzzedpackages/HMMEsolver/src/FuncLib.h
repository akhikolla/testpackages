
#ifndef _HMMEsolver_FUNCLIB_H
#define _HMMEsolver_FUNCLIB_H


//#include <iostream>
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
using namespace arma;


mat GetWtilde(double* W, double* lambda, int n, int m);

mat GetXWs(mat X, double* _W, double* _s);
mat GetZWs(double* W, double* s, double* psi, double* lambda, int n, int m);
mat GetZWX(mat X, double* _W, int n, int m);


bool CppGetDiag(double* out, double* X, double* Y, double* mu, double* lambda, int* n, int* m, int*p);

#endif


