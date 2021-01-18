/*
 * cblas_zhpr.c
 * The program is a C interface to zhpr.
 *
 * Keita Teranishi  3/23/98
 *
 */

#ifndef CBLAS_ZHPR_H_
#define CBLAS_ZHPR_H_

#include <stdio.h>
#include <stdlib.h>

inline void cblas_zhpr(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                       const int N, const double alpha, const void *X,
                       const int incX, void *A) {
  char UL;
#ifdef F77_CHAR
  F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif

#ifdef F77_INT
  F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incx
#endif
  // int n, i, tincx, incx=incX;
  int n, i, tincx;
  double *x = (double *)X, *xx = (double *)X, *tx, *st;

  extern int CBLAS_CallFromC;
  extern int RowMajorStrg;
  RowMajorStrg = 0;

  CBLAS_CallFromC = 1;
  if (layout == CblasColMajor) {
    if (Uplo == CblasLower)
      UL = 'L';
    else if (Uplo == CblasUpper)
      UL = 'U';
    else {
      cblas_xerbla(2, "cblas_zhpr", "Illegal Uplo setting, %d\n", Uplo);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }
#ifdef F77_CHAR
    F77_UL = C2F_CHAR(&UL);
#endif

    F77_NAME(zhpr)
    (F77_UL, &F77_N, &alpha, (const Rcomplex *)X, &F77_incX, (Rcomplex *)A);

  } else if (layout == CblasRowMajor) {
    RowMajorStrg = 1;
    if (Uplo == CblasUpper)
      UL = 'L';
    else if (Uplo == CblasLower)
      UL = 'U';
    else {
      cblas_xerbla(2, "cblas_zhpr", "Illegal Uplo setting, %d\n", Uplo);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }
#ifdef F77_CHAR
    F77_UL = C2F_CHAR(&UL);
#endif
    if (N > 0) {
      n = N << 1;
      x = (double *)malloc(n * sizeof(double));
      tx = x;
      if (incX > 0) {
        i = incX << 1;
        tincx = 2;
        st = x + n;
      } else {
        i = incX * (-2);
        tincx = -2;
        st = x - 2;
        x += (n - 2);
      }
      do {
        *x = *xx;
        x[1] = -xx[1];
        x += tincx;
        xx += i;
      } while (x != st);
      x = tx;
#ifdef F77_INT
      F77_incX = 1;
#else
      incx = 1;
#endif
    } else
      x = (double *)X;

    F77_NAME(zhpr)
    (F77_UL, &F77_N, &alpha, (const Rcomplex *)x, &F77_incX, (Rcomplex *)A);

  } else {
    cblas_xerbla(1, "cblas_zhpr", "Illegal layout setting, %d\n", layout);
    CBLAS_CallFromC = 0;
    RowMajorStrg = 0;
    return;
  }
  if (X != x) free(x);
  CBLAS_CallFromC = 0;
  RowMajorStrg = 0;
  return;
}

#endif
