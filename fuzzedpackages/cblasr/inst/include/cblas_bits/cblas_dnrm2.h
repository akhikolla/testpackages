/*
 * cblas_dnrm2.c
 *
 * The program is a C interface to dnrm2.
 * It calls the fortranwrapper before calling dnrm2.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DNRM2_H_
#define CBLAS_DNRM2_H_

inline double cblas_dnrm2(const int N, const double *X, const int incX) {
  double nrm2;
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  nrm2 = F77_NAME(dnrm2)(&F77_N, X, &F77_incX);
  return nrm2;
}

#endif
