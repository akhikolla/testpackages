/*
 * cblas_dscal.c
 *
 * The program is a C interface to dscal.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DSCAL_H_
#define CBLAS_DSCAL_H_

inline void cblas_dscal(const int N, const double alpha, double *X,
                        const int incX) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  F77_NAME(dscal)(&F77_N, &alpha, X, &F77_incX);
}

#endif
