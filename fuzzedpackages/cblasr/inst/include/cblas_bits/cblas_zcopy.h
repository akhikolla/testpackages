/*
 * cblas_zcopy.c
 *
 * The program is a C interface to zcopy.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZCOPY_H_
#define CBLAS_ZCOPY_H_

inline void cblas_zcopy(const int N, const void *X, const int incX, void *Y,
                        const int incY) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
  F77_NAME(zcopy)
  (&F77_N, (const Rcomplex *)X, &F77_incX, (const Rcomplex *)Y, &F77_incY);
}

#endif
