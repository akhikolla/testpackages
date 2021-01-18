
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>  
#define NODEBUG

/* Useful constants                        */ 
const static double SQR3 = 1.73205080756888;
const static double SQR5 = 2.23606797749979;
const static double EPSILON = 1e-7;

/*============================================================================*  
 * AUTHORS                                                                    *
 *                                                                            *
 * Yves Deville, Olivier Roustant and David Ginsbourger for the 'ReDice'      *
 * consortium www.redice-project.org                                          *
 *                                                                            *
 * DESCRIPTION                                                                *
 *                                                                            *
 * Kernels with one-dimensional site args intended to be used through         *
 * the .Call interface.                                                       *
 *                                                                            *
 * o Several couples of sites can be computed in a loop: 'x1' and 'x2' point  *
 *   to vectors having the same length 'n'.                                   *
 *                                                                            *
 * o The gradient of the kernel (with respect to the parameters in 'par') is  *
 *   returned as a matrix  attribute of the kernel value with one column by   *
 *   parameter derivation.                                                    *
 *                                                                            *
 * NOTE                                                                       *
 *                                                                            *
 * The variance is given as the last element of the 'par' vector.             *
 * Hence the derivation with respect to this parameter is formally done (the  *
 * derivative is the correlation kernel).                                     *
 *                                                                            *
 *============================================================================*/

/*============================================================================*
 * exp. Note that the derivatives for x = 0 does not exists but is set to     *
 * zero, which is suitable for gradient computations.                         *
 *============================================================================*/

SEXP k1FunExpC(SEXP x) {   /* x = site                                       */  

  int i, n = LENGTH(x);
  SEXP value, der, attrNm, dd;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(value = allocVector(REALSXP, n));  
  PROTECT(dd = allocVector(INTSXP, 2));
  INTEGER(dd)[0] = n;
  INTEGER(dd)[1] = 2; 
  PROTECT(der = allocArray(REALSXP, dd));      

  double z, emz, *rx = REAL(x), *rvalue = REAL(value), *rder = REAL(der);

  for (i = 0; i < n; i++) {
    z = fabs(rx[i]);
    emz = exp(-z);
    rvalue[i] = emz;
    if (z < EPSILON) {
      rder[i] = 0.0;
      rder[n + i] = 0.0;
    } else {
      rder[i] = -emz;
      rder[n + i] = emz;
    }
  }
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("der"));
  SET_ATTR(value, attrNm, der);
  
  UNPROTECT(5);
  return(value);

}


/*============================================================================*
 * Matern 3/2                                                                 *
 *============================================================================*/

SEXP k1FunMatern3_2C(SEXP x) {   /* x = site                                       */  

  int i, n = LENGTH(x);
  SEXP value, der, attrNm, dd;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(value = allocVector(REALSXP, n));  
  PROTECT(dd = allocVector(INTSXP, 2));
  INTEGER(dd)[0] = n;
  INTEGER(dd)[1] = 2; 
  PROTECT(der = allocArray(REALSXP, dd));      

  double z, emz, *rx = REAL(x), *rvalue = REAL(value), *rder = REAL(der);

  for (i = 0; i < n; i++) {
    z = SQR3 * fabs(rx[i]);
    emz = exp(-z);
    rvalue[i] = (1.0 + z) * emz;
    rder[i] = -SQR3 * z* emz;
    rder[n + i] = 3.0 * (-1.0 + z) * emz;
  }
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("der"));
  SET_ATTR(value, attrNm, der);
  
  UNPROTECT(5);
  return(value);

}

/*============================================================================*
 * Matern 3/2                                                                 *
 *============================================================================*/

SEXP k1FunMatern5_2C(SEXP x) {   /* x = site                                       */  

  int i, n = LENGTH(x);
  SEXP value, der, attrNm, dd;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(value = allocVector(REALSXP, n));  
  PROTECT(dd = allocVector(INTSXP, 2));
  INTEGER(dd)[0] = n;
  INTEGER(dd)[1] = 2; 
  PROTECT(der = allocArray(REALSXP, dd));      

  double z, emz, *rx = REAL(x), *rvalue = REAL(value), *rder = REAL(der);

  for (i = 0; i < n; i++) {
    z = SQR5 * fabs(rx[i]);
    emz = exp(-z);
    rvalue[i] = (1.0 + z * ( 1.0 + z / 3.0) ) * emz;
    rder[i] = -SQR5 * z * (1.0 + z) *  emz / 3.0;
    rder[n + i] = -SQR5 * (1.0 + z * (1.0 - z)) * emz / 3.0;
  }
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("der"));
  SET_ATTR(value, attrNm, der);
  
  UNPROTECT(5);
  return(value);

}

/*============================================================================*
 * Gauss                                                                      *
 *============================================================================*/

SEXP k1FunGaussC(SEXP x) {   /* x = site                                       */  

  int i, n = LENGTH(x);
  SEXP value, der, attrNm, dd;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(value = allocVector(REALSXP, n));  
  PROTECT(dd = allocVector(INTSXP, 2));
  INTEGER(dd)[0] = n;
  INTEGER(dd)[1] = 2; 
  PROTECT(der = allocArray(REALSXP, dd));      

  double *rx = REAL(x), *rvalue = REAL(value), *rder = REAL(der);

  for (i = 0; i < n; i++) {
    rvalue[i] = exp(-rx[i] * rx[i] / 2.0);
    rder[i] = -rx[i] * rvalue[i];
    rder[n + i] = (-1.0 + rx[i] * rx[i]) * rvalue[i];
  }
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("der"));
  SET_ATTR(value, attrNm, der);
  
  UNPROTECT(5);
  return(value);

}

/*============================================================================*
 * Power Exponential. Caution: with a parameter!                              *
 *============================================================================*/

SEXP k1FunPowExpC(SEXP x, SEXP alpha) {   /* x = site                          */  

  int i, n = LENGTH(x);
  SEXP value, der, grad, attrNm, dd;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(value = allocVector(REALSXP, n));  
  PROTECT(dd = allocVector(INTSXP, 2));
  INTEGER(dd)[0] = n;
  INTEGER(dd)[1] = 2; 
  PROTECT(der = allocArray(REALSXP, dd));
  INTEGER(dd)[1] = 1; 
  PROTECT(grad = allocArray(REALSXP, dd));
  
  double z, *rx = REAL(x), *rvalue = REAL(value),
    *rder = REAL(der), *rgrad = REAL(grad), *ralpha = REAL(alpha);
  
  for (i = 0; i < n; i++) {
    z = pow(fabs(rx[i]), ralpha[0]); 
    rvalue[i] = exp(-z);
    if (z > EPSILON) {    
      rder[i] = -ralpha[0] * z * rvalue[i] / rx[i];
      rder[n + i] = - (ralpha[0] * (z - 1.0) + 1.0) * rder[i] / rx[i];
      rgrad[i] = -log(fabs(rx[i])) * z * rvalue[i];
    }
    else {
      rder[i] = 0.0;
      rder[n + i] = 0.0;
      rgrad[i] = 0.0;
    }
  }
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("der"));
  SET_ATTR(value, attrNm, der);
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  SET_ATTR(value, attrNm, grad);

  UNPROTECT(7);
  return(value);

}

