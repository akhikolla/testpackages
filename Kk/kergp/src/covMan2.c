#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#define NODEBUG

/*============================================================================*
 * AUTHOR                                                                     *
 *                                                                            *
 *    Yves Deville <deville.yves@alpestat.com> for the ReDice Consortium      * 
 * www.redice-project.org.                                                    *
 *                                                                            *
 * DESCRIPTION                                                                *
 *                                                                            *
 * Compute the contributions 'gamma_k' of the derivarives w.r.t. the parame-  *
 * -ter r in the log-Likelihood. More precisely, this is                      *
 *                                                                            *
 *       gamma_k = sum_{i,j} W_{ij} \partial_r C_{ij}                         *
 *                                                                            *
 * where W is symetric matrix with dimension n * n. Only the lower part of    *
 * W is passed here as a vector 'weights' of length n * (n + 1) /2. Diagonal  *
 * elements must have formerly been  be divided by 2 in order to get the 1/2  *
 * factor in the log-likelihood.                                              *
 *                                                                            *
 *                                                                            *
 *============================================================================*/

SEXP scores_covMan(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		    SEXP Xt,       // Transpose of the spatial design : d*n
		    SEXP par,      // vector of 'npar' param for the Kd kernel
		    SEXP weights,  // Vector of length m := n * (n + 1) / 2
		    SEXP rho) {    // An R environment

  int  i, j, ij, k, n, m, ipar, d, npar;
  
  double *rxt = REAL(Xt), *rx1, *rx2,
    *rweights = REAL(weights), *rscores;
  
  SEXP dimXt, R_fcall, x1, x2,
    kernValue, dkernValue,  attrNm, scores;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(Xt)) error("'Xt' must be a matrix");
  if(!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  PROTECT(Xt = coerceVector(Xt, REALSXP));
  PROTECT(dimXt = getAttrib(Xt, R_DimSymbol));
  d = INTEGER(dimXt)[0];
  n = INTEGER(dimXt)[1]; 
  m = ( n * (n + 1) ) / 2;
  
#ifdef DEBUG 
  Rprintf("d = %d n = %d, m = %d\n", d, n, m);
#endif 
  
  /* check the parameters arrays  */
  PROTECT(par = coerceVector(par, REALSXP));
  npar = LENGTH(par);

#ifdef DEBUG
  double *rpar = REAL(par)
  for (ipar = 0; ipar < npar; ipar++) {
    Rprintf("par[%d] = %7.3f  ", ipar, rpar[ipar]);
  }
  Rprintf("\n");
#endif

  PROTECT(weights = coerceVector(weights, REALSXP));
  if (LENGTH(weights) != m) {
    error("vector 'weights' on input with bad length");
  }
  
  PROTECT(x1 = allocVector(REALSXP, d)); //  sites x1, x2
  PROTECT(x2 = allocVector(REALSXP, d));
  PROTECT(scores = allocVector(REALSXP, npar));

  rx1 = REAL(x1);
  rx2 = REAL(x2);
  rscores = REAL(scores);
  
  PROTECT(R_fcall = lang4(fun, x1, x2, par));
  SETCADDDR(R_fcall, par);

  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

  // initialise the scores. 'ipar' is an index in the whole param vector
  for (ipar = 0; ipar < npar; ipar++) {
    rscores[ipar] = 0.0;
  }
 
  ij = 0;
  
  for (i = 0; i < n; i++) { 
    // copy
    for (k = 0; k < d; k++) {
      rx1[k] = rxt[i*d + k];
    }
    SETCADR(R_fcall, x1);
    
    for (j = i; j < n; j++) {
      // copy
      for (k = 0; k < d; k++) {
	rx2[k] = rxt[j*d + k]; 
      }
      SETCADDR(R_fcall, x2);
      kernValue = eval(R_fcall, rho);
      // Rprintf("Kernel value = %7.3f\n", REAL(kernValue)[0]);
      dkernValue = GET_ATTR(kernValue, attrNm);
      // for (ipar = 0; ipar < npar; ipar++) {
      // Rprintf("%d dK = %7.3f\n", ipar, REAL(dkernValue)[ipar]);
      //}
      //Rprintf("Kernel value = %7.3f\n", REAL(kernValue)[1]);
      // update the scores
      for (ipar = 0; ipar < npar; ipar++) {
	rscores[ipar] += rweights[ij] * REAL(dkernValue)[ipar];
      }
      ij++;
    }
  }

  UNPROTECT(9);
  return(scores);

}
