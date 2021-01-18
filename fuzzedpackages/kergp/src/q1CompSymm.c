#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>  
#define NODEBUG 

/*============================================================================*
  
  GOAL

  Compute a correlation matrix or its Cholesky lower root using the
  compound symmetry pattern.

  AUTHOR 

  Yves Deville <deville.yves@alpestat.com>

  DETAILS
  
  We do not consider here the case where 'lowerSQRT' is FALSE, because
  it is straightforward in R, and no gain would result from using the
  '.Call' API in this case.

 *============================================================================*/

SEXP corLev_CompSymm(SEXP par, 
		     SEXP nlevels,
		     SEXP lowerSQRT,
		     SEXP compGrad) {
  
  int i, j, j1, npar = LENGTH(par), m = INTEGER(nlevels)[0],  
    m1 = m - 1, m2 = m * m;

  if (npar != 1) error("length of 'par' not equal to 1");
  if (!INTEGER(lowerSQRT)[0])  error("'lowerSQRT must be TRUE");
  
  SEXP cor;
  
  PROTECT(par = coerceVector(par, REALSXP));
  double *rpar = REAL(par);

  PROTECT(cor = allocMatrix(REALSXP, m, m));
  double *rcor = REAL(cor);

  /*=======================================================================
    a[j] is the diagonal element in column j, b[j] is the common
    value of the elements below the diagonal in column j.
    ======================================================================= */
  
  double *a = (double *) R_alloc(m, sizeof(double));
  double *b = (double *) R_alloc(m1, sizeof(double));
  a[0] = 1.0;
  rcor[0] = 1.0;
  double S2 = 0.0;

  if (INTEGER(compGrad)[0]) {

    SEXP dcor, attrNm;

    PROTECT(dcor = allocVector(REALSXP, m2));
    PROTECT(attrNm = NEW_CHARACTER(1));

    double *rdcor = REAL(dcor);
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
     
    /* Initialise */ 
    
    for (i = 0; i < m; i++)  {
      for (j = 0; j < m; j++)  {
	rcor[i + j * m] = 0.0;
	rdcor[i + j * m] = 0.0;
      }
    }
    
    rcor[0] = 1.0;

    /*=======================================================================
      Compute the derivative of the elements of 'a' and 'b' w.r.t. 'rho'
      by recursion.
      ======================================================================= */

    double *da = (double *) R_alloc(m, sizeof(double));
    double *db = (double *) R_alloc(m1, sizeof(double));
    double dS2 = 0.0;
    da[0] = 0.0;

    for (j = 0; j < m1; j++) {
      
      b[j] = (rpar[0] - S2) / a[j];
      db[j] = (1.0 - dS2 - b[j] * da[j]) / a[j];

      S2 += b[j] * b[j];
      dS2 += 2 * b[j] * db[j];

      j1 = j + 1;
      a[j1] = sqrt(1.0 - S2);
      rcor[j1 + m * j1] = a[j1];

      da[j1] = -dS2 / a[j1] / 2.0;
      rdcor[j1 + m * j1] = da[j1];

      for (i = j1; i < m; i++) {
	rcor[i + j * m] = b[j];
	rdcor[i + j * m] = db[j];
      }

    }

    SET_ATTR(cor, attrNm, dcor);
    UNPROTECT(4);
    return(cor);

  } else {

     /* Initialise */ 
    
    for (i = 0; i < m; i++)  {
      for (j = 0; j < m; j++)  {
	rcor[i + j * m] = 0.0;
      }
    }

    rcor[0] = 1.0;

    for (j = 0; j < m1; j++) {
      
      b[j] = (rpar[0] - S2) / a[j];
      S2 += b[j] * b[j];

      j1 = j + 1;
      a[j1] = sqrt(1.0 - S2);
      rcor[j1 + m * j1] = a[j1];

      for (i = j1; i < m; i++) {
	rcor[i + j * m] = b[j];
      }
      
    }

    UNPROTECT(2);  
    return(cor);
  
  }

}
