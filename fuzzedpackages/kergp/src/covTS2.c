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
 * elements must have formerly been  be divided by 2 in order to get the 1/2 
 * factor in the log-likelihood.                                              *
 *                                                                            *
 * CAUTION                                                                    *
 * -------                                                                    *
 * At the time, it is implicitely assumed that for each dimension 'ell', at   *
 * most one parameter in 'par_ell' depends on the derivation parameter in     *
 * 'par'. Thus a PowExp kernel can not have is shape and its scale set to     *
 * the same parameter. This assumption is easy to drop.                       * 
 *                                                                            *
 *============================================================================*/

SEXP scores_covTS(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		  SEXP Xt,       // Transpose of the spatial design : d*n
		  SEXP par,      // vector of 'npar' param for the Kd kernel
		  SEXP parMap,   // TRANSPOSE of the Kd mapping : d * npar
		  SEXP weights,  // Vector of length m := n * (n + 1) / 2
		  SEXP rho) {    // An R environment
  
  int  i, j, ij, k, ell, n, m, d, p, npar, *iparMap = INTEGER(parMap),  
    ipoint;
 
  double *rxt = REAL(Xt), *rx1_ell, *rx2_ell, *S,
    *rpar = REAL(par), *rpar_ell, *rweights = REAL(weights), *rscores;
  
  SEXP dimXt, dimParMap, R_fcall, x1_ell, x2_ell, par_ell,
    kernValue, dkernValue,  attrNm, scores;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(Xt)) error("'Xt' must be a matrix");
  if(!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  Xt = coerceVector(Xt, REALSXP);
  PROTECT(dimXt = getAttrib(Xt, R_DimSymbol));
  d = INTEGER(dimXt)[0];
  n = INTEGER(dimXt)[1]; 
  m = ( n * (n + 1) ) / 2;
  
#ifdef DEBUG 
  Rprintf("d = %d n = %d, m = %d\n", d, n, m);
#endif 
  
  /* check the parameters arrays  */
  npar = LENGTH(par);
  par = coerceVector(par, REALSXP);
  parMap = coerceVector(parMap, INTSXP);  // length npar * d
  PROTECT(dimParMap = getAttrib(parMap, R_DimSymbol));
  p = INTEGER(dimParMap)[0];

  S = (double *) R_alloc(npar, sizeof(double));

#ifdef DEBUG
  Rprintf("Total number of parameters npar = %d\n", npar);
  Rprintf("Number of the one-dimensional kernel parameters p = %d\n", p);
  Rprintf("Mapping of the parameters\n");
  for (ell = 0; ell < d; ell++) {
    for (k = 0; k < p; k++) {
      Rprintf("   ell = %d, k = %d, parMap[ell, k] = %d\n", 
	      ell, k, iparMap[ell * p + k]);
    }
  }
#endif 
  
  weights = coerceVector(weights, REALSXP);
  if (LENGTH(weights) != m) {
    error("vector 'weights' on input with bad length");
  }
  
  PROTECT(x1_ell = allocVector(REALSXP, 1)); // one-dim sites x1, x2
  PROTECT(x2_ell = allocVector(REALSXP, 1));
  PROTECT(par_ell = allocVector(REALSXP, p));
  PROTECT(scores = allocVector(REALSXP, npar));

  rx1_ell = REAL(x1_ell);
  rx2_ell = REAL(x2_ell);
  rpar_ell = REAL(par_ell);
  rscores = REAL(scores);

  PROTECT(R_fcall = lang4(fun, x1_ell, x2_ell, par_ell)); 
  PROTECT(kernValue = allocVector(REALSXP, 1));
  PROTECT(dkernValue = allocVector(REALSXP, p));
  
  PROTECT(attrNm = NEW_CHARACTER(1));
  SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

  // initialise the scores. 'ipoint' is an index in the whole param vector
  for (ipoint = 0; ipoint < npar; ipoint++) {
    rscores[ipoint] = 0.0;
  }
 
  ij = 0;
  
  for (i = 0; i < n; i++) { 
    for (j = i; j < n; j++) { 
      // initialise S. 'ipoint' is an index in the whole param vector.
      for (ipoint = 0; ipoint < npar; ipoint++) {
	S[ipoint] = 0.0;
      }
      // loop over spatail dimensions
      for (ell = 0; ell < d; ell++) {
	
	rx1_ell[0] = rxt[i * d + ell];
	SETCADR(R_fcall, x1_ell);
	
	rx2_ell[0] = rxt[j * d + ell];
	SETCADDR(R_fcall, x2_ell);
	
	for (k = 0; k < p; k++) {
	  ipoint = iparMap[ell * p + k];
	  rpar_ell[k] = rpar[ipoint];
	}
	
	SETCADDDR(R_fcall, par_ell);
	kernValue = eval(R_fcall, rho);
	dkernValue = GET_ATTR(kernValue, attrNm);
	
	for (k = 0; k < p; k++) {
	  ipoint = iparMap[ell * p + k];
	  S[ipoint] += REAL(dkernValue)[k];
	}
	
      }

      // for (ipoint = 0; ipoint < npar; ipoint++) {   
      // Rprintf("S[%d] = %4.3f\n", ipoint, S[ipoint]);
      // }
      
      for (ipoint = 0; ipoint < npar; ipoint++) {
	rscores[ipoint] += rweights[ij] * S[ipoint];
      }

      ij++;

    }
  }

  UNPROTECT(10);
  return(scores);
 
}
