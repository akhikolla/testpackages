/*
 * cblas_izamax.c
 *
 * The program is a C interface to izamax.
 * It calls the fortran wrapper before calling izamax.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_IZAMAX_H_
#define CBLAS_IZAMAX_H_

inline CBLAS_INDEX cblas_izamax(const int N, const void *X, const int incX) {
  CBLAS_INDEX iamax;
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  iamax = F77_NAME(izamax)(&F77_N, (const Rcomplex *)X, &F77_incX);
  return (iamax ? iamax - 1 : 0);
}

#endif
