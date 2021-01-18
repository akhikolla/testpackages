/*
 * cblas_zaxpy.c
 *
 * The program is a C interface to zaxpy.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZAXPY_H_
#define CBLAS_ZAXPY_H_

inline void cblas_zaxpy(const int N, const void *alpha, const void *X,
                        const int incX, void *Y, const int incY) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
  F77_NAME(zaxpy)
  (&F77_N, (const Rcomplex *)alpha, (const Rcomplex *)X, &F77_incX,
   (const Rcomplex *)Y, &F77_incY);
}

#endif
