/*
 * cblas_idamax.c
 *
 * The program is a C interface to idamax.
 * It calls the fortran wrapper before calling idamax.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_IDAMAX_H_
#define CBLAS_IDAMAX_H_

inline CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX) {
  CBLAS_INDEX iamax;
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  iamax = F77_NAME(idamax)(&F77_N, X, &F77_incX);
  return iamax ? iamax - 1 : 0;
}

#endif
