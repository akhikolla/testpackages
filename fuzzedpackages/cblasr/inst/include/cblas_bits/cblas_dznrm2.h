/*
 * cblas_dznrm2.c
 *
 * The program is a C interface to dznrm2.
 * It calls the fortran wrapper before calling dznrm2.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DZNRM2_H_
#define CBLAS_DZNRM2_H_

inline double cblas_dznrm2(const int N, const void *X, const int incX) {
  double nrm2;
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  nrm2 = F77_NAME(dznrm2)(&F77_N, (const Rcomplex *)X, &F77_incX);
  return nrm2;
}

#endif
