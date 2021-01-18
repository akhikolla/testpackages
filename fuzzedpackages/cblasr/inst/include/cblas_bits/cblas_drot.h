/*
 * cblas_drot.c
 *
 * The program is a C interface to drot.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DROT_H_
#define CBLAS_DROT_H_

inline void cblas_drot(const int N, double *X, const int incX, double *Y,
                       const int incY, const double c, const double s) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
  F77_NAME(drot)(&F77_N, X, &F77_incX, Y, &F77_incY, &c, &s);
  return;
}

#endif
