#ifndef _HMMEsolver_MAINLIB_H
#define _HMMEsolver_MAINLIB_H

#ifdef __cplusplus
extern "C" {
#endif

  void SolveHMM(double* out, double* X, double* Y, double* mu, double* lambda, int* n, int* m, int* p);


#ifdef __cplusplus
}
#endif


#endif
