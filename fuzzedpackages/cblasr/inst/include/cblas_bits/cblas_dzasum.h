/*
 * cblas_dzasum.c
 *
 * The program is a C interface to dzasum.
 * It calls the fortran wrapper before calling dzasum.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DZASUM_H_
#define CBLAS_DZASUM_H_

inline double cblas_dzasum(const int N, const void *X, const int incX) {
  double asum;
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  asum = F77_NAME(dzasum)(&F77_N, (const Rcomplex *)X, &F77_incX);
  return asum;
}

#endif
