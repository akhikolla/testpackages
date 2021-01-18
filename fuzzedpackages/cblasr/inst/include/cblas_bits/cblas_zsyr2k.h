/*
 *
 * cblas_zsyr2k.c
 * This program is a C interface to zsyr2k.
 * Written by Keita Teranishi
 * 4/8/1998
 *
 */

#ifndef CBLAS_ZSYR2K_H_
#define CBLAS_ZSYR2K_H_

inline void cblas_zsyr2k(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                         const CBLAS_TRANSPOSE Trans, const int N, const int K,
                         const void *alpha, const void *A, const int lda,
                         const void *B, const int ldb, const void *beta,
                         void *C, const int ldc) {
  char UL, TR;
#ifdef F77_CHAR
  F77_CHAR F77_TR, F77_UL;
#else
#define F77_TR &TR
#define F77_UL &UL
#endif

#ifdef F77_INT
  F77_INT F77_N = N, F77_K = K, F77_lda = lda, F77_ldb = ldb;
  F77_INT F77_ldc = ldc;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
#endif

  extern int CBLAS_CallFromC;
  extern int RowMajorStrg;
  RowMajorStrg = 0;
  CBLAS_CallFromC = 1;

  if (layout == CblasColMajor) {
    if (Uplo == CblasUpper)
      UL = 'U';
    else if (Uplo == CblasLower)
      UL = 'L';
    else {
      cblas_xerbla(2, "cblas_zsyr2k", "Illegal Uplo setting, %d\n", Uplo);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }

    if (Trans == CblasTrans)
      TR = 'T';
    else if (Trans == CblasConjTrans)
      TR = 'C';
    else if (Trans == CblasNoTrans)
      TR = 'N';
    else {
      cblas_xerbla(3, "cblas_zsyr2k", "Illegal Trans setting, %d\n", Trans);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR(&UL);
    F77_TR = C2F_CHAR(&TR);
#endif

    F77_NAME(zsyr2k)
    (F77_UL, F77_TR, &F77_N, &F77_K, (Rcomplex *)alpha, (Rcomplex *)A, &F77_lda,
     (Rcomplex *)B, &F77_ldb, (Rcomplex *)beta, (Rcomplex *)C, &F77_ldc);
  } else if (layout == CblasRowMajor) {
    RowMajorStrg = 1;
    if (Uplo == CblasUpper)
      UL = 'L';
    else if (Uplo == CblasLower)
      UL = 'U';
    else {
      cblas_xerbla(3, "cblas_zsyr2k", "Illegal Uplo setting, %d\n", Uplo);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }
    if (Trans == CblasTrans)
      TR = 'N';
    else if (Trans == CblasConjTrans)
      TR = 'N';
    else if (Trans == CblasNoTrans)
      TR = 'T';
    else {
      cblas_xerbla(3, "cblas_zsyr2k", "Illegal Trans setting, %d\n", Trans);
      CBLAS_CallFromC = 0;
      RowMajorStrg = 0;
      return;
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR(&UL);
    F77_TR = C2F_CHAR(&TR);
#endif

    F77_NAME(zsyr2k)
    (F77_UL, F77_TR, &F77_N, &F77_K, (Rcomplex *)alpha, (Rcomplex *)A, &F77_lda,
     (Rcomplex *)B, &F77_ldb, (Rcomplex *)beta, (Rcomplex *)C, &F77_ldc);
  } else
    cblas_xerbla(1, "cblas_zsyr2k", "Illegal layout setting, %d\n", layout);
  CBLAS_CallFromC = 0;
  RowMajorStrg = 0;
  return;
}

#endif
