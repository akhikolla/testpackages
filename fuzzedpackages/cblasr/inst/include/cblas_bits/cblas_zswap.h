/*
 * cblas_zswap.c
 *
 * The program is a C interface to zswap.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZSWAP_H_
#define CBLAS_ZSWAP_H_

inline void cblas_zswap(const int N, void *X, const int incX, void *Y,
                        const int incY) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
  F77_NAME(zswap)(&F77_N, (Rcomplex *)X, &F77_incX, (Rcomplex *)Y, &F77_incY);
}

#endif
