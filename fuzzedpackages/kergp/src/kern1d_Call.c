
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
 * Exponential                                                                *
 *============================================================================*/

SEXP k1ExpC(SEXP x1,       /* site #1 (scalar)                                 */
	   SEXP x2,       /* site #2 (scalar)                                 */
	   SEXP par) {    /* parameter vector: c(range, var)                  */  

  int i1, i2, n1 = LENGTH(x1), n2 = LENGTH(x2);
  SEXP kernValue, dkernValue, attrNm, dd;
  
  PROTECT(x1 = coerceVector(x1, REALSXP));
  PROTECT(x2 = coerceVector(x2, REALSXP));
  PROTECT(par = coerceVector(par, REALSXP));
  
  if (LENGTH(par) != 2) {
    Rprintf("length(par) = %d\n", LENGTH(par)); 
    error("For \"Exp\" kernel, 'par' must be of length 2");
  }

  double z, kStar,  *rk, *rdk,
    *rx1 = REAL(x1), *rx2 = REAL(x2), *rpar = REAL(par);

  PROTECT(dd = allocVector(INTSXP, 3));
  PROTECT(kernValue = allocMatrix(REALSXP, n1, n2));        
  INTEGER(dd)[0] = n1; INTEGER(dd)[1] = n2; INTEGER(dd)[2] = 2; 
  PROTECT(dkernValue = allocArray(REALSXP, dd)); 
  rk = REAL(kernValue);
  rdk =  REAL(dkernValue);

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      z = fabs(rx1[i1] - rx2[i2]) / rpar[0];
      kStar = exp(-z);
      rk[i1 + i2 * n1] = kStar * rpar[1];
      rdk[i1 + i2 * n1 ] = z * kStar * rpar[1] / rpar[0];
      rdk[i1 + i2 * n1 + n1 * n2] = kStar;
    }
  }

  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  // attribut "gradient" de 'Cov'.
  SET_ATTR(kernValue, attrNm, dkernValue);
  
  UNPROTECT(7);
  return(kernValue);

}

/*============================================================================*
 * Gauss                                                                      *
 *============================================================================*/

SEXP k1GaussC(SEXP x1,       /* site #1 (scalar)                               */
	     SEXP x2,       /* site #2 (scalar)                               */
	     SEXP par) {    /* parameter vector: c(range, var)                */  

  int i1, i2, n1 = LENGTH(x1), n2 = LENGTH(x2);
  SEXP kernValue, dkernValue, attrNm, dd;
  
  PROTECT(x1 = coerceVector(x1, REALSXP));
  PROTECT(x2 = coerceVector(x2, REALSXP));
  PROTECT(par = coerceVector(par, REALSXP));

  if (LENGTH(par) != 2) {
    error("For \"Gauss\" kernel, 'par' must be of length 2");
  }

  double z, z2, kStar, *rk, *rdk,
    *rx1 = REAL(x1), *rx2 = REAL(x2), *rpar = REAL(par);

  PROTECT(dd = allocVector(INTSXP, 3));
  PROTECT(kernValue = allocMatrix(REALSXP, n1, n2));    
  INTEGER(dd)[0] = n1; INTEGER(dd)[1] = n2; INTEGER(dd)[2] = 2; 
  PROTECT(dkernValue = allocArray(REALSXP, dd));      

  rk = REAL(kernValue);
  rdk =  REAL(dkernValue);

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      z = fabs(rx1[i1] - rx2[i2]) / rpar[0];
      z2 = z * z;
      kStar = exp(-z2 / 2.0);
      rk[i1 + i2 * n1] = kStar * rpar[1];
      rdk[i1 + i2 * n1] = z2 * kStar * rpar[1] / rpar[0]; 
      rdk[i1 + i2 * n1 + n1 * n2] = kStar;
    }
  }

  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  // attribut "gradient" de 'Cov'.
  SET_ATTR(kernValue, attrNm, dkernValue);
  
  UNPROTECT(7);
  return(kernValue);

}

/*============================================================================*
 *  Power-exponential                                                         *
 *  Note the check z != 0.0 for the derivative                                *
 *============================================================================*/

SEXP k1PowExpC(SEXP x1,       /* site #1 (scalar)                              */
	      SEXP x2,       /* site #2 (scalar)                              */
	      SEXP par) {    /* parameter vector: c(range, shape, var)        */  
  
  int i1, i2, n1 = LENGTH(x1), n2 = LENGTH(x2);
  SEXP kernValue, dkernValue, attrNm, dd;
  
  PROTECT(x1 = coerceVector(x1, REALSXP));
  PROTECT(x2 = coerceVector(x2, REALSXP));
  PROTECT(par = coerceVector(par, REALSXP));

  if (LENGTH(par) != 3) {
    error("For \"PowExp\" kernel, 'par' must be of length 3");
  }

  double z, zAlpha, kStar, *rk, *rdk,
    *rx1 = REAL(x1), *rx2 = REAL(x2), *rpar = REAL(par);

  PROTECT(dd = allocVector(INTSXP, 3));
  PROTECT(kernValue = allocMatrix(REALSXP, n1, n2));    
  INTEGER(dd)[0] = n1; INTEGER(dd)[1] = n2; INTEGER(dd)[2] = 3;      
  PROTECT(dkernValue = allocArray(REALSXP, dd)); 
  rk = REAL(kernValue);
  rdk =  REAL(dkernValue);

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      z = fabs(rx1[i1] - rx2[i2]) / rpar[0];
      zAlpha = pow(z, rpar[1]); 
      kStar = exp(-zAlpha);
      rk[i1 + i2 * n1] = kStar * rpar[2];
      rdk[i1 + i2 * n1] =  zAlpha * kStar * rpar[1] * rpar[2] / rpar[0];
      if (z < EPSILON) {
	rdk[i1 + i2 * n1 + n1 * n2] = 0.0;
      }  
      else {
	rdk[i1 + i2 * n1 + n1 * n2] = -log(z) * zAlpha * kStar * rpar[2];
      }
      rdk[i1 + i2 * n1 + 2 * n1 * n2] = kStar;
    }
  }

  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  // attribut "gradient" de 'Cov'.
  SET_ATTR(kernValue, attrNm, dkernValue);
  
  UNPROTECT(7);
  return(kernValue);

}

/*============================================================================*
 *  Matern nu = 3/2                                                           *
 *============================================================================*/

SEXP k1Matern3_2C(SEXP x1,       /* site #1 (scalar)                           */
		 SEXP x2,       /* site #2 (scalar)                           */
		 SEXP par) {    /* parameter vector: c(range, var)            */  
   
  int i1, i2, n1 = LENGTH(x1), n2 = LENGTH(x2);
  SEXP kernValue, dkernValue, attrNm, dd;
  
  PROTECT(x1 = coerceVector(x1, REALSXP));
  PROTECT(x2 = coerceVector(x2, REALSXP));
  PROTECT(par = coerceVector(par, REALSXP));
  
  if (LENGTH(par) != 2) {
    error("For \"Matern3_2\" kernel, 'par' must be of length 2");
  }
  
  double z, kStar, expz, *rk, *rdk, 
    *rx1 = REAL(x1), *rx2 = REAL(x2), *rpar = REAL(par);
   
  PROTECT(dd = allocVector(INTSXP, 3));
  PROTECT(kernValue = allocMatrix(REALSXP, n1, n2));    
  INTEGER(dd)[0] = n1; INTEGER(dd)[1] = n2; INTEGER(dd)[2] = 2;      
  PROTECT(dkernValue = allocArray(REALSXP, dd)); 

  rk = REAL(kernValue);
  rdk =  REAL(dkernValue);
  
  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {

      z = fabs(rx1[i1] - rx2[i2]) * SQR3 / rpar[0];  // z = sqrt(3) h / theta
      expz = exp(-z);
      kStar = expz * (1.0 + z);
      rk[i1 + i2 * n1] = kStar * rpar[1];
      rdk[i1 + i2 * n1] = z * z * expz * rpar[1] / rpar[0] ; 
      rdk[i1 + i2 * n1 + n1 * n2] = kStar;
    }
  }
    
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  // attribut "gradient" de 'Cov'.
  SET_ATTR(kernValue, attrNm, dkernValue);
  
  UNPROTECT(7);
  return(kernValue);


} 

/*============================================================================*
 *  Matern nu = 5/2                                                           *
 *============================================================================*/

SEXP k1Matern5_2C(SEXP x1,       /* site #1 (scalar)                           */
		 SEXP x2,       /* site #2 (scalar)                           */
		 SEXP par) {    /* parameter vector: c(range,  var)           */
  
  int i1, i2, n1 = LENGTH(x1), n2 = LENGTH(x2);
  SEXP kernValue, dkernValue, attrNm, dd;
  
  PROTECT(x1 = coerceVector(x1, REALSXP));
  PROTECT(x2 = coerceVector(x2, REALSXP));
  PROTECT(par = coerceVector(par, REALSXP));

  if (LENGTH(par) != 2) {
    error("For \"Matern5_2\" kernel, 'par' must be of length 2");
  }
  
  double z, kStar, expz, *rk, *rdk,
    *rx1 = REAL(x1), *rx2 = REAL(x2), *rpar = REAL(par);
   
  PROTECT(dd = allocVector(INTSXP, 3));
  PROTECT(kernValue = allocMatrix(REALSXP, n1, n2));    
  INTEGER(dd)[0] = n1; INTEGER(dd)[1] = n2; INTEGER(dd)[2] = 2;      
  PROTECT(dkernValue = allocArray(REALSXP, dd));

  rk = REAL(kernValue);
  rdk =  REAL(dkernValue);
  
  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      z = fabs(rx1[i1] - rx2[i2]) * SQR5 / rpar[0];
      expz = exp(-z);
      kStar = (1.0 + z * (1.0 + z / 3.0)) * expz; 
      rk[i1 + i2 * n1] = kStar * rpar[1];
      rdk[i1 + i2 * n1] = z * z * (1.0 + z) * expz * rpar[1] / rpar[0] / 3.0; 
      rdk[i1 + i2 * n1 + n1 * n2] = kStar;
    }
  }

  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
  // attribut "gradient" de 'Cov'.
  SET_ATTR(kernValue, attrNm, dkernValue);
  
  UNPROTECT(7);
  return(kernValue);

}
