/*
 * cblas_zdotu_sub.c
 *
 * The program is a C interface to zdotu.
 * It calls the fortran wrapper before calling zdotu.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_ZDOTU_SUB_H_
#define CBLAS_ZDOTU_SUB_H_

inline void cblas_zdotu_sub(const int N, const void *X, const int incX,
                            const void *Y, const int incY, void *dotu) {
#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
  *((Rcomplex *)dotu) = F77_NAME(zdotu)(&F77_N, (const Rcomplex *)X, &F77_incX,
                                        (const Rcomplex *)Y, &F77_incY);
  return;
}

#endif
