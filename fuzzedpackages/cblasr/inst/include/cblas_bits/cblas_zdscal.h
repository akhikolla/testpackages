/*
 * cblas_zdscal.c
 *
 * The program is a C interface to zdscal.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZDSCAL_H_
#define CBLAS_ZDSCAL_H_

inline void cblas_zdscal(const int N, const double alpha, void *X,
                         const int incX) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
  F77_NAME(zdscal)(&F77_N, &alpha, (Rcomplex *)X, &F77_incX);
}

#endif
