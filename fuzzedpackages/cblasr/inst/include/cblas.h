#ifndef CBLAS_H
#define CBLAS_H

#include <stddef.h>

#include <R.h>
#include <R_ext/BLAS.h>

#ifdef __cplusplus
extern "C" { /* Assume C declarations for C++ */
#endif       /* __cplusplus */

/*
 * Enumerated and derived types
 */
#ifdef WeirdNEC
#define CBLAS_INDEX long
#else
#define CBLAS_INDEX int
#endif

typedef enum { CblasRowMajor = 101, CblasColMajor = 102 } CBLAS_LAYOUT;
typedef enum {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
} CBLAS_TRANSPOSE;
typedef enum { CblasUpper = 121, CblasLower = 122 } CBLAS_UPLO;
typedef enum { CblasNonUnit = 131, CblasUnit = 132 } CBLAS_DIAG;
typedef enum { CblasLeft = 141, CblasRight = 142 } CBLAS_SIDE;

typedef CBLAS_LAYOUT
    CBLAS_ORDER; /* this for backward compatibility with CBLAS_ORDER */

#define F77_INT int

#include "cblas_bits/cblas_globals.h"
#include "cblas_bits/cblas_xerbla.h"
#include "cblas_bits/xerbla.h"

/* Double Precision Level 1 BLAS */

#include "cblas_bits/cblas_dasum.h"
#include "cblas_bits/cblas_daxpy.h"
#include "cblas_bits/cblas_dcopy.h"
#include "cblas_bits/cblas_ddot.h"
#include "cblas_bits/cblas_dnrm2.h"
#include "cblas_bits/cblas_drot.h"
#include "cblas_bits/cblas_drotg.h"
#include "cblas_bits/cblas_drotm.h"
#include "cblas_bits/cblas_drotmg.h"
#include "cblas_bits/cblas_dscal.h"
#include "cblas_bits/cblas_dswap.h"
#include "cblas_bits/cblas_idamax.h"

/* Double Precision Level 2 BLAS */

#include "cblas_bits/cblas_dgbmv.h"
#include "cblas_bits/cblas_dgemv.h"
#include "cblas_bits/cblas_dsbmv.h"
#include "cblas_bits/cblas_dspmv.h"
#include "cblas_bits/cblas_dsymv.h"
#include "cblas_bits/cblas_dtbmv.h"
#include "cblas_bits/cblas_dtpmv.h"
#include "cblas_bits/cblas_dtrmv.h"

#include "cblas_bits/cblas_dtbsv.h"
#include "cblas_bits/cblas_dtpsv.h"
#include "cblas_bits/cblas_dtrsv.h"

#include "cblas_bits/cblas_dger.h"
#include "cblas_bits/cblas_dspr.h"
#include "cblas_bits/cblas_dspr2.h"
#include "cblas_bits/cblas_dsyr.h"
#include "cblas_bits/cblas_dsyr2.h"

/* Double Precision Level 3 BLAS */

#include "cblas_bits/cblas_dgemm.h"
#include "cblas_bits/cblas_dtrmm.h"
#include "cblas_bits/cblas_dtrsm.h"

#include "cblas_bits/cblas_dsymm.h"
#include "cblas_bits/cblas_dsyr2k.h"
#include "cblas_bits/cblas_dsyrk.h"

/* Double complex BLAS routines added for 2.3.0 */

// dcabs1?
#include "cblas_bits/cblas_dzasum.h"
#include "cblas_bits/cblas_dznrm2.h"
#include "cblas_bits/cblas_izamax.h"
#include "cblas_bits/cblas_zaxpy.h"
#include "cblas_bits/cblas_zcopy.h"

/*
 * WARNING from R_ext/BLAS.h! The next two return a value that may not
 * be compatible between C and Fortran, and even if it is, this might
 * not be the right translation to C.  Only use after
 * configure-testing with your compilers.
 */
#include "cblas_bits/cblas_zdotc_sub.h"
#include "cblas_bits/cblas_zdotu_sub.h"

// zdrot?
#include "cblas_bits/cblas_zdscal.h"
#include "cblas_bits/cblas_zgbmv.h"
#include "cblas_bits/cblas_zgemm.h"
#include "cblas_bits/cblas_zgemv.h"
#include "cblas_bits/cblas_zgerc.h"
#include "cblas_bits/cblas_zgeru.h"
#include "cblas_bits/cblas_zhbmv.h"
#include "cblas_bits/cblas_zhemm.h"
#include "cblas_bits/cblas_zhemv.h"
#include "cblas_bits/cblas_zher.h"
#include "cblas_bits/cblas_zher2.h"
#include "cblas_bits/cblas_zher2k.h"
#include "cblas_bits/cblas_zherk.h"
#include "cblas_bits/cblas_zhpmv.h"
#include "cblas_bits/cblas_zhpr.h"
#include "cblas_bits/cblas_zhpr2.h"
// zrotg?
#include "cblas_bits/cblas_zscal.h"
#include "cblas_bits/cblas_zswap.h"
#include "cblas_bits/cblas_zsymm.h"
#include "cblas_bits/cblas_zsyr2k.h"
#include "cblas_bits/cblas_zsyrk.h"
#include "cblas_bits/cblas_ztbmv.h"
#include "cblas_bits/cblas_ztbsv.h"
#include "cblas_bits/cblas_ztpmv.h"
#include "cblas_bits/cblas_ztpsv.h"
#include "cblas_bits/cblas_ztrmm.h"
#include "cblas_bits/cblas_ztrmv.h"
#include "cblas_bits/cblas_ztrsm.h"
#include "cblas_bits/cblas_ztrsv.h"

#ifdef __cplusplus
}
#endif

#endif
