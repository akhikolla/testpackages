#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>  
#define NODEBUG 

/*============================================================================*
   
  NOTE

  This program is part of the 'kergp' R package.

  GOAL

  Compute a correlation matrix or its Cholesky lower root using the
  general "spherical" (or "Cholesky") parameterisation of Pinheiro and
  Bates.

  AUTHOR 

  Yves Deville <deville.yves@alpestat.com> adapting existing codes by
  Douglas M. Bates Jose C. Pinheiro (nlme package), and by Eperan
  Padonou.
  
  CAUTION 

  The parameters are the angles between O and pi. THEY ARE ASSUMED TO
  BE GIVEN IN ROW ORDER theta[1, 1], theta[2, 1], theta[2, 2],
  theta[3, 1], ...

  DETAILS
      
  'ell' is an index referring to the 'par' array, so 'ell' ranges from
  0 to m * (m - 1) / 2.  'ell_k' and ell_i are other such indices used
  to cope with elements in row 'i'.

 *============================================================================*/

SEXP corLev_Symm(SEXP par, 
		 SEXP nlevels,
		 SEXP lowerSQRT,
		 SEXP compGrad) {

  int i, j, k, npar = LENGTH(par), m = INTEGER(nlevels)[0], 
    m2 = m * m, ell, ell_i, ell_k, ell_m =  m * (m - 1) / 2;

  if (npar != ell_m) {
    error("length of 'par' not equal to 'm * (m - 1) / 2'");
  }
  
  double aux, eps = 1e-8;
  
  /* 'Dbg' is used in place of compiler directive to avoid
     indentation bugs with emacs. */
  
  int Dbg = 0;
  SEXP cor;

  PROTECT(par = coerceVector(par, REALSXP));
  double *rpar = REAL(par);

  PROTECT(cor = allocMatrix(REALSXP, m, m));
  double *rcor = REAL(cor);

  if (Dbg) {  
    Rprintf("compGrad = %d, lowerSQRT = %d\n", INTEGER(compGrad)[0], 
	    INTEGER(lowerSQRT)[0]);
    Rprintf("m = %d, npar = %d, m2 = %d\n", m, npar, m2);
  }

  if (INTEGER(compGrad)[0]) {

    int *ells = (int *) R_alloc(m + 1, sizeof(int));
    ells[m] =  ell_m;
    double *c = (double *) R_alloc(m, sizeof(double));
    double *s = (double *) R_alloc(m, sizeof(double));
 
    SEXP dcor, attrNm;
    
    PROTECT(dcor = allocVector(REALSXP, m2 * npar));
    PROTECT(attrNm = NEW_CHARACTER(1));
    double *rdcor = REAL(dcor);
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

   /*========================================================================
     First compute the lower Cholesky root of the covariance matrix
    =========================================================================*/

    /* Initialise */ 
    
    if (Dbg)  Rprintf("Initialise\n");
    
    for (i = 0; i < m; i++)  {
      for (j = 0; j < m; j++)  {
	rcor[i + j * m] = 0.0;
	for (ell = 0; ell < npar; ell++) {
	  rdcor[i + j * m + ell * m2] = 0.0;
	}
      }
    }
    
    /* loop on rows */ 
    ell = 0;
    
    for (i = 0; i < m; i++)  {
      
      ell_i = ell;
      ells[i] = ell;
      aux = 1.0;
      
      if (Dbg) Rprintf("row i = %d ell_i = %d\n", i, ell_i);
      
      for (j = 0; j < i; j++) {
	
	c[j] = cos(rpar[ell]);
	s[j] = sin(rpar[ell]);
	
	rcor[i + j * m] =  c[j] * aux; 
	aux *= s[j]; 
	
	/* derivatives for k < j  */
	ell_k = ell_i;
	
	for (k = 0; k < j; k++) { 
	  
	  if (s[k] > eps) { 
	    rdcor[i + j * m + ell_k * m2] =  rcor[i + j * m] * c[k] / s[k];
	    if (Dbg) {
	      Rprintf("   deriv. i = %d, j = %d w.r.t. par[%d]\n", i, j, ell_k);
	    }
	  } 
	  
	  ell_k++;
	  
	}
	
	/* now derivative for k = j.  Note that the following could
	   work when the angle theta_{ij} is not pi / 2.  

	   rdcor[i + j * m + ell_k * m2] = -rcor[i + j * m] * s[j] / c[j];  

	   but is less direct.
	*/
	rdcor[i + j * m + ell_k * m2] = - aux;
	if (Dbg) {
	  Rprintf("   deriv. i = %d, j = %d w.r.t. par[%d]\n", i, j, ell_k);
	}
	
	ell++;
	
      }
      
      /* diagonal element */
      rcor[i + i * m] =  aux;     
      
      /* derivatives for k < i  */
      ell_k = ell_i;
      
      for (k = 0; k < i; k++) { 	
	
	if (s[k] > eps) { 
	  rdcor[i + i * m + ell_k * m2] =  rcor[i + i * m] * c[k] / s[k];
	  if (Dbg) {
	    Rprintf("   deriv. i = %d, j = %d w.r.t. par[%d]\n", i, j, ell_k);
	  }
	} 
	ell_k++;
	
      }
      
    }

    /*========================================================================
      Now compute the correlation matrix from its lower Cholesky square root
      ========================================================================*/
    
    if (!INTEGER(lowerSQRT)[0]) {
      
      double *temp = (double *) R_alloc(m2, sizeof(double));
      double *dtemp = (double *) R_alloc(m2 * npar, sizeof(double));
      
      if (Dbg) Rprintf("Build the correlations and their derivatives\n");
      
      for (i = 0; i < m; i++) {
	for (j = 0; j < m; j++) {
	  for (ell = 0; ell < npar; ell++) {
	    dtemp[i + j * m + ell * m2] = 0.0;
	  } 
	  temp[i + j * m] = 0.0;
	}
      }

      /*==================================================================
	Compute the product  LL'
	==================================================================*/
      
      for (i = 0; i < m; i++) {
	for (j = 0; j <= i; j++) {
	  aux = 0.0;
	  for (k = 0; k <= j; k++) {
	    aux += rcor[i + k * m] * rcor[j + k * m]; 
	  }
	  temp[i + j * m] = aux;
	  temp[j + i * m] = aux;
	}
      }
	
      for (i = 0; i < m; i++) {
	for (j = 0; j <= i; j++) {	  
	  for (k = 0; k <= j; k++) {
	    
	    /*================================================================
	      the following two loops do what a single huge loop for
	      ell = 0 to npar would do, but remind that most of
	      derivatives in 'rdcor' are zero.
	      ================================================================*/
	    
	    ell_k = ells[j] + k;
	    if (ell_k >= ell_m) ell_k = ell_m - 1;

	    for (ell = ells[j]; ell <= ell_k; ell++) {
	      dtemp[i + j * m + ell * m2] +=  rcor[i + k * m] * 
		rdcor[j + k * m + ell * m2];
	    }

	    ell_k = ells[i] + k;
	    if (ell_k >= ell_m) ell_k = ell_m - 1;

	    for (ell = ells[i]; ell <= ell_k; ell++) {
	      dtemp[i + j * m + ell * m2] +=  rcor[j + k * m] * 
		rdcor[i + k * m + ell * m2];
	    }
	    
	  }
	}
      }

      /*================================================================
	Copy. XXX could be made faster : copy only non-zero elements!
	However, copying is not as expensive as is computing.
       ==================================================================*/

      for (i = 0; i < m; i++) {
	for (j = 0; j < m; j++) {
	  rcor[i + j * m] = temp[i + j * m];
	  for (ell = 0; ell < npar; ell++) {
	    rdcor[i + j * m + ell * m2] = dtemp[i + j * m + ell * m2]; 
	    rdcor[j + i * m + ell * m2] = dtemp[i + j * m + ell * m2];
	  }
	}
      }

    }

    if (Dbg) Rprintf("Done: exit\n");

    SET_ATTR(cor, attrNm, dcor);
    UNPROTECT(4);
    
    return(cor);
      
  } else {
    
    /*=======================================================================
      Now consider the much easier no-gradient case
      =======================================================================*/
    
    /* Initialise */ 
    if (Dbg) Rprintf("Initialise\n");
    
    for (i = 0; i < m; i++)  {
      for (j = 0; j < m; j++)  {
	rcor[i + j * m] = 0.0;
      }
    }

    /* loop on rows */ 
    ell = 0;
    
    for (i = 0; i < m; i++)  {
      
      aux = 1.0;
      
      for (j = 0; j < i; j++) {	
	rcor[i + j * m] =  cos(rpar[ell]) * aux; 
	aux *= sin(rpar[ell]); 
	  ell++;
      }
      
      rcor[i + i * m] = aux;
      
    }
    
    if (!INTEGER(lowerSQRT)[0]) {
      
      /* Build the correlations. Could be done using a BLAS routine */ 

      double *temp = (double *) R_alloc(m2, sizeof(double));
      
      if (Dbg) Rprintf("Build correlations\n");
      
      for (i = 0; i < m; i++) {
	for (j = 0; j <= i; j++) {
	  aux = 0.0;
	  for (k = 0; k <= j; k++) {
	    aux += rcor[i + k * m] * rcor[j + k * m]; 
	  }
	  temp[i + j * m] = aux;
	  temp[j + i * m] = aux;
	}
      }
      
      /* copy */
      for (i = 0; i < m; i++) {
	for (j = 0; j < m; j++) {
	  rcor[i + j * m] = temp[i + j * m]; 
	}
      }
      
    }
    
    UNPROTECT(2);
    return(cor);
    
  }
   
}
