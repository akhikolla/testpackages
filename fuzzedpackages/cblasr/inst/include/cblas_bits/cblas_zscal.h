/*
 * cblas_zscal.c
 *
 * The program is a C interface to zscal.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZSCAL_H_
#define CBLAS_ZSCAL_H_

inline void cblas_zscal(const int N, const void *alpha, void *X,
                        const int incX) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  F77_NAME(zscal)(&F77_N, (const Rcomplex *)alpha, (Rcomplex *)X, &F77_incX);
}

#endif
