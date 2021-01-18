
#include "FuncLib.h"
#include "MainLib.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>


void SolveHMM(double* out, double* X, double* Y, double* mu, double* lambda, int* n, int* m, int* p) {

  bool Result;
  Result = CppGetDiag(out, X, Y, mu, lambda, n, m, p);

  return;

}

