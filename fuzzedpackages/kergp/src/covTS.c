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
 * Compute (/update in the future)  the covariance matrix 'Cov' by adding the *
 * contribution of a 'covTS'  additive term.                                  * 
 *                                                                            *
 *                                                                            *
 *       sum_{ell = 1}^d  K_1( x_{i, ell}, x_{j, ell}; theta_ell )            *
 *                                                                            *
 * where K1 is some one-dimensional covariance kernel depending on a vector   *
 * 'theta' of 'p' parameters. The total number of parameter for the d-dimen-  *
 * -sional kernel is 'npar' and falls between 1 and  p * d.                   *
 *                                                                            *
 * In this version, the argument 'par' can be used and passed to the kernel   *
 * function from the C-level. Moreover, the gradient of the matrix with       *
 * respect to ONE of the parameters can be computed if wanted, and then       *
 * returned as the "gradient" attribute of the covariance. Note that          *
 * "gradient is NOT the full gradient (which is an huge array) but only one   *
 * 2-dimensional  slice of it. The argument 'index' is the index of the para- *
 * -meter with respect to which the derivation is done. It is in C-fashion      *
 * (index >=0) and differs from the argument of the R function by one unit.   *
 *                                                                            *
 * The covariance should be passed as an argument for additive models, and    *
 * initialized to zeroes from the R-code if necessary.                        *
 *                                                                            *
 * NOTE                                                                       *
 *                                                                            *
 * At the time, it is implicitely assumed that for each dimension 'ell', at   *
 * most one parameter in 'par_ell' depends on the derivation parameter in     *
 * 'par'. Thus a PowExp kernel can not have is shape and its scale set to     *
 * the same parameter. This assumption is easy to drop.                       * 
 *                                                                            *
 *============================================================================*/

SEXP covMat_covTS(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		  SEXP Xt,       // Transpose of the spatial design : d*n
		  SEXP par,      // vector of 'npar' param for the Kd kernel
		  SEXP parMap,   // TRANSPOSE of the Kd mapping : d * npar
		  SEXP compGrad, // Integer 0 or 1: compute gradient matrix?
		  SEXP index,    // Index for grad.
		  SEXP rho) {    // An R environment
  
  int  i, j, k, ell, n, d, p, *iparMap = INTEGER(parMap),  ipoint;
  
  double *rxt = REAL(Xt), *rx1_ell, *rx2_ell, *rCov, 
    *rpar = REAL(par), *rpar_ell;
  
  SEXP dimXt, dimParMap, R_fcall, Cov, x1_ell, x2_ell, par_ell;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(Xt)) error("'Xt' must be a matrix");
  if (!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  PROTECT(Xt = coerceVector(Xt, REALSXP));
  PROTECT(dimXt = getAttrib(Xt, R_DimSymbol));
  d = INTEGER(dimXt)[0];
  n = INTEGER(dimXt)[1]; 

#ifdef DEBUG 
  Rprintf("d = %d n = %d\n", d, n);
#endif

  /* check the parameters arrays  */
  par = coerceVector(par, REALSXP);
  parMap = coerceVector(parMap, INTSXP);  // length npar * d
  PROTECT(dimParMap = getAttrib(parMap, R_DimSymbol));
  p = INTEGER(dimParMap)[0];

#ifdef DEBUG
  int npar = LENGTH(par);
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

  // PROTECT(Cov = allocMatrix(REALSXP, n, n)); //allocate the n x n matrix 
  PROTECT(Cov = allocMatrix(REALSXP, n, n));
  PROTECT(x1_ell = allocVector(REALSXP, 1)); // one-dim sites x1, x2
  PROTECT(x2_ell = allocVector(REALSXP, 1));
  PROTECT(par_ell = allocVector(REALSXP, p));

  rCov = REAL(Cov);
  rx1_ell = REAL(x1_ell);
  rx2_ell = REAL(x2_ell);
  rpar_ell = REAL(par_ell);

  PROTECT(R_fcall = lang4(fun, x1_ell, x2_ell, par_ell));

  // SETCADDDR(R_fcall, par);
  // SETCAD4R(R_fcall, compGrad);    

  // FOR FUTURE IMPLEMENTATION: grad flag could be passed to the 1d kernel in
  // order to save time when no derivation is required.
  
  // initialize upper triangle j >= i
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) { 
      rCov[i + j * n] = 0.0;
    }
  }
 
  if ( INTEGER(compGrad)[0] ) {

    SEXP dCov, kernValue, dkernValue, attrNm;
    double *rdCov;
    int iindex, ipointGrad;

    iindex = INTEGER(index)[0];

#ifdef DEBUG 
    Rprintf("Gradient required for index %d\n", iindex);
#endif

    PROTECT(dCov = allocMatrix(REALSXP, n, n)); //allocate the n x n matrix
    PROTECT(kernValue = allocVector(REALSXP, 1));
    PROTECT(dkernValue = allocVector(REALSXP, p));
    
    PROTECT(attrNm = NEW_CHARACTER(1));
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    rdCov = REAL(dCov);
    
    // initialize upper triangle j >= i
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) { 
	rdCov[i + j * n] = 0.0;
      }
    }

    for (ell = 0; ell < d; ell++) {

#ifdef DEBUG 
      Rprintf("\nDimension ell = %d\n", ell);
#endif

      // o select parameter vector for dim 'ell',
      // o watch for a required derivation (ipoint == iindex))
      // o fill the vector 'rpar_ell' representing 'theta_ell'
      ipointGrad = -1;

      for (k = 0; k < p; k++) {
	ipoint = iparMap[ell * p + k];
	if (ipoint == iindex) ipointGrad = k;
	rpar_ell[k] = rpar[ipoint];
      }

#ifdef DEBUG 
      if (ipointGrad >= 0) {
	Rprintf("   derivation param used as index: %d \n", ipointGrad);
      } else {
	Rprintf("    no deriv. for this input %d\n", ipointGrad);
      }
#endif

      SETCADDDR(R_fcall, par_ell);

      if (ipointGrad >= 0) {
	// the parameter for dim 'ell' depends on par[iindex]: derivation.
	for (i = 0; i < n; i++) {
	  rx1_ell[0] = rxt[i * d + ell];
	  SETCADR(R_fcall, x1_ell);

	  for (j = i; j < n; j++) {
	    rx2_ell[0] = rxt[j * d + ell];
	    SETCADDR(R_fcall, x2_ell);
	    kernValue = eval(R_fcall, rho);
	    rCov[i + j * n] += REAL(kernValue)[0];
	    dkernValue = GET_ATTR(kernValue, attrNm);
	    rdCov[i + j * n] += REAL(dkernValue)[ipointGrad];
	    //  Rprintf("  grad = %8.6f\n", REAL(dkernValue)[ipointGrad]);
	  }

	}

      } else {
	for (i = 0; i < n; i++) {
	  // the parameter for dim 'ell' does not depend on  par[iindex]
	  rx1_ell[0] = rxt[i * d + ell];
	  SETCADR(R_fcall, x1_ell);

	  for (j = i; j < n; j++) {
	    rx2_ell[0] = rxt[j * d + ell];
	    SETCADDR(R_fcall, x2_ell);
	    kernValue = eval(R_fcall, rho);
	    rCov[i + j * n] += REAL(kernValue)[0];    
	  }
	}
      }
    }
#ifdef DEBUG   
    Rprintf("Copy in lower triangle and exit\n\n");
#endif
    // copy both lower triangles j < i
    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++) { 
	rCov[i + j * n] = rCov[j + i * n];
	rdCov[i + j * n] = rdCov[j + i * n];
      }
    }

    // set the gradient attribute of 'Cov'.
    SET_ATTR(Cov, attrNm, dCov);
    UNPROTECT(12);
    return(Cov);
    
  } else {
#ifdef DEBUG   
    Rprintf("No gradient required\n");
#endif

    SEXP kernValue;
    PROTECT(kernValue = allocVector(REALSXP, 1));

    for (ell = 0; ell < d; ell++) {
#ifdef DEBUG  
      Rprintf("Dimension ell = %d\n", ell);
#endif
      // select parameter vector for dim 'ell' 
      for (k = 0; k < p; k++) {
	ipoint = iparMap[ell * p + k];
	rpar_ell[k] = rpar[ipoint];
#ifdef DEBUG 
	Rprintf("k = %d rpar_ell[k] = %7.4f\n", k, rpar_ell[k]);
#endif
      }
      SETCADDDR(R_fcall, par_ell);
      for (i = 0; i < n; i++) {
	rx1_ell[0] = rxt[i * d + ell];
	SETCADR(R_fcall, x1_ell);
	for (j = i; j < n; j++) {
	  rx2_ell[0] = rxt[j * d + ell];
	  SETCADDR(R_fcall, x2_ell);
	  kernValue = eval(R_fcall, rho);
	  rCov[i + j * n] += REAL(kernValue)[0];
	}
      }
      
    }
#ifdef DEBUG 
    Rprintf("Copy in lower triangle and exit \n\n");
#endif
    // copy lower triangle j < i
    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++) { 
	rCov[i + j * n] = rCov[j + i * n];
      }
    }
    UNPROTECT(9);
    return(Cov);
    
  }

}

/*============================================================================*
 * covMat1Mat2                                                                *
 *                                                                            *
 * At the time, no derivation is possible.                                    * 
 *                                                                            *
 *============================================================================*/

SEXP covMatMat_covTS(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		     SEXP X1t,      // Transpose of the spatial design : d*n1
		     SEXP X2t,      // Transpose of the spatial design : d*n2
		     SEXP par,      // vector of 'npar' param for the Kd kernel
		     SEXP parMap,   // TRANSPOSE of the Kd mapping : d * npar
		     SEXP compGrad, // NOT USED
		     SEXP index,    // NOT USED
		     SEXP rho) {    // An R environment
  
  int  i, j, k, ell, n1, n2, d, p, *iparMap = INTEGER(parMap), ipoint;
  
  double *rx1t = REAL(X1t),  *rx2t = REAL(X2t),
    *rx1_ell, *rx2_ell, *rCov, 
    *rpar = REAL(par), *rpar_ell;
  
  SEXP dimX1t, dimX2t, dimParMap, R_fcall, Cov, x1_ell, x2_ell, par_ell;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(X1t)) error("'X1t' must be a matrix"); 
  if (!isMatrix(X2t)) error("'X2t' must be a matrix");
  if (!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  PROTECT(X1t = coerceVector(X1t, REALSXP));
  PROTECT(dimX1t = getAttrib(X1t, R_DimSymbol));
  d = INTEGER(dimX1t)[0];
  n1 = INTEGER(dimX1t)[1];
 
  PROTECT(X2t = coerceVector(X2t, REALSXP));
  PROTECT(dimX2t = getAttrib(X2t, R_DimSymbol));
  if (INTEGER(dimX2t)[0] != d) {
    error("'X1t' and 'X2t must have the same number of rows (number of inputs)");
  }
  n2 = INTEGER(dimX2t)[1]; 

#ifdef DEBUG 
  Rprintf("d = %d n1 = %d n2 =%d\n", d, n1);
#endif

  /* check the parameters arrays  */
  par = coerceVector(par, REALSXP);
  parMap = coerceVector(parMap, INTSXP);  // length npar * d
  PROTECT(dimParMap = getAttrib(parMap, R_DimSymbol));
  p = INTEGER(dimParMap)[0];

#ifdef DEBUG
  int npar = LENGTH(par);
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

  // PROTECT(Cov = allocMatrix(REALSXP, n, n)); //allocate the n x n matrix 
  PROTECT(Cov = allocMatrix(REALSXP, n1, n2));
  PROTECT(x1_ell = allocVector(REALSXP, 1)); // one-dim sites x1, x2
  PROTECT(x2_ell = allocVector(REALSXP, 1));
  PROTECT(par_ell = allocVector(REALSXP, p));

  rCov = REAL(Cov);
  rx1_ell = REAL(x1_ell);
  rx2_ell = REAL(x2_ell);
  rpar_ell = REAL(par_ell);

  PROTECT(R_fcall = lang4(fun, x1_ell, x2_ell, par_ell));

  // FOR FUTURE IMPLEMENTATION: grad flag could be passed to the 1d kernel in
  // order to save time when no derivation is required.
  
  // initialize at zero
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) { 
      rCov[i + j * n1] = 0.0;
    }
  }
 
  if ( INTEGER(compGrad)[0] ) {

    error("Gradient computation not implemented for covMatMat");
    
  } else {

    SEXP kernValue;
    PROTECT(kernValue = allocVector(REALSXP, 1));

    for (ell = 0; ell < d; ell++) {

#ifdef DEBUG  
      Rprintf("Dimension ell = %d\n", ell);
#endif
      // select parameter vector for dim 'ell' 
      for (k = 0; k < p; k++) {
	ipoint = iparMap[ell * p + k];
	rpar_ell[k] = rpar[ipoint];
#ifdef DEBUG 
	Rprintf("k = %d rpar_ell[k] = %7.4f\n", k, rpar_ell[k]);
#endif
      }
      SETCADDDR(R_fcall, par_ell);
      for (i = 0; i < n1; i++) {
	rx1_ell[0] = rx1t[i * d + ell];
	SETCADR(R_fcall, x1_ell);
	for (j = 0; j < n2; j++) {
	  rx2_ell[0] = rx2t[j * d + ell];
	  SETCADDR(R_fcall, x2_ell);
	  kernValue = eval(R_fcall, rho);
	  rCov[i + j * n1] += REAL(kernValue)[0];
	}
      }
      
    }

    UNPROTECT(11);
    return(Cov);
    
  }

}

/*============================================================================*
 * AUTHOR                                                                     *
 *                                                                            *
 *    Yves Deville <deville.yves@alpestat.com>                                *
 *                                                                            *
 * DESCRIPTION                                                                *
 *                                                                            *
 * Compute (/update in the future)  the variance vector 'Var' by adding the   *
 * contribution of a 'covTS'  additive term.                                  * 
 *                                                                            *
 *                                                                            *
 *       sum_{ell = 1}^d  K_1( x_{i, ell}, x_{i, ell}; theta_ell )            *
 *                                                                            *
 * where K1 is some one-dimensional covariance kernel depending on a vector   *
 * 'theta' of 'p' parameters. The total number of parameter for the d-dimen-  *
 * -sional kernel is 'npar' and falls between 1 and  p * d.                   *
 *                                                                            *
 * See the explanations for the function 'covMat_covTS' in this file          *
 *============================================================================*/

SEXP varVec_covTS(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		  SEXP Xt,       // Transpose of the spatial design : d*n
		  SEXP par,      // vector of 'npar' param for the Kd kernel
		  SEXP parMap,   // TRANSPOSE of the Kd mapping : d * npar
		  SEXP compGrad, // Integer 0 or 1: compute gradient matrix?
		  SEXP index,    // Index for grad.
		  SEXP rho) {    // An R environment
  
  int  i, k, ell, n, d, p, *iparMap = INTEGER(parMap),  ipoint;
  
  double *rxt = REAL(Xt), *rx1_ell, *rVar, 
    *rpar = REAL(par), *rpar_ell;
  
  SEXP dimXt, dimParMap, R_fcall, Var, x1_ell, par_ell;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(Xt)) error("'Xt' must be a matrix");
  if (!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  PROTECT(Xt = coerceVector(Xt, REALSXP));
  PROTECT(dimXt = getAttrib(Xt, R_DimSymbol));
  d = INTEGER(dimXt)[0];
  n = INTEGER(dimXt)[1]; 

  /* check the parameters arrays  */
  par = coerceVector(par, REALSXP);
  parMap = coerceVector(parMap, INTSXP);  // length npar * d
  PROTECT(dimParMap = getAttrib(parMap, R_DimSymbol));
  p = INTEGER(dimParMap)[0];

  PROTECT(Var = allocVector(REALSXP, n));
  PROTECT(x1_ell = allocVector(REALSXP, 1)); // one-dim sites x1
  PROTECT(par_ell = allocVector(REALSXP, p));

  rVar = REAL(Var);
  rx1_ell = REAL(x1_ell);
  rpar_ell = REAL(par_ell);

  PROTECT(R_fcall = lang4(fun, x1_ell, x1_ell, par_ell));
  
  // FOR FUTURE IMPLEMENTATION: grad flag could be passed to the 1d kernel in
  // order to save time when no derivation is required.
  
  // initialize 
  for (i = 0; i < n; i++) { 
    rVar[i] = 0.0;
  }
  
  if ( INTEGER(compGrad)[0] ) {
    
    SEXP dVar, kernValue, dkernValue, attrNm;
    double *rdVar;
    int iindex, ipointGrad;

    iindex = INTEGER(index)[0];
    
    PROTECT(dVar = allocMatrix(REALSXP, n, 1)); //allocate the n x 1 matrix
    PROTECT(kernValue = allocVector(REALSXP, 1));
    PROTECT(dkernValue = allocVector(REALSXP, p));
    
    PROTECT(attrNm = NEW_CHARACTER(1));
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    rdVar = REAL(dVar);
    
    // initialize
    for (i = 0; i < n; i++) { 
	rdVar[i] = 0.0;
    }

    for (ell = 0; ell < d; ell++) {

      // o select parameter vector for dim 'ell',
      // o watch for a required derivation (ipoint == iindex))
      // o fill the vector 'rpar_ell' representing 'theta_ell'
      ipointGrad = -1;

      for (k = 0; k < p; k++) {
	ipoint = iparMap[ell * p + k];
	if (ipoint == iindex) ipointGrad = k;
	rpar_ell[k] = rpar[ipoint];
      }

      SETCADDDR(R_fcall, par_ell);

      if (ipointGrad >= 0) {
	// the parameter for dim 'ell' depends on par[iindex]: derivation.
	for (i = 0; i < n; i++) {
	  rx1_ell[0] = rxt[i * d + ell];
	  SETCADR(R_fcall, x1_ell);
	  SETCADDR(R_fcall, x1_ell);
	  kernValue = eval(R_fcall, rho);
	  rVar[i] += REAL(kernValue)[0];
	  dkernValue = GET_ATTR(kernValue, attrNm);
	  rdVar[i] += REAL(dkernValue)[ipointGrad];

	}

      } else {
	for (i = 0; i < n; i++) {
	  // the parameter for dim 'ell' does not depend on  par[iindex]
	  rx1_ell[0] = rxt[i * d + ell];
	  SETCADR(R_fcall, x1_ell);
	  SETCADDR(R_fcall, x1_ell);
	  kernValue = eval(R_fcall, rho);
	  rVar[i] += REAL(kernValue)[0];    
	  
	}
      }
    }
    // set the gradient attribute of 'Var'.
    SET_ATTR(Var, attrNm, dVar);
    UNPROTECT(11);
    return(Var);
    
  } else {

    for (ell = 0; ell < d; ell++) {
      // select parameter vector for dim 'ell' 
      for (k = 0; k < p; k++) {
	ipoint = iparMap[ell * p + k];
	rpar_ell[k] = rpar[ipoint];
      }
      SETCADDDR(R_fcall, par_ell);
      for (i = 0; i < n; i++) {
	rx1_ell[0] = rxt[i * d + ell];
	SETCADR(R_fcall, x1_ell);
	SETCADDR(R_fcall, x1_ell);
	rVar[i] += REAL(eval(R_fcall, rho))[0];
      }
      
    }
    UNPROTECT(7);
    return(Var);
    
  }

}
