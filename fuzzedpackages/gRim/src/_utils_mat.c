/* ************************************************************
   
   Various matrix utility functions
   
   Author: Søren Højsgaard
   
   ************************************************************ */ 

#define USE_FC_LEN_T
#include <Rconfig.h>

#include <Rdefines.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif
//#include <R_ext/Applic.h>       /* for dgemm */
//#include <R_ext/RS.h>		/* for F77_... */

#include "_utils_print.h"

#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP))


/* ************************************************************
   ************************************************************
   ************************************************************
   DECLARATIONS
   ************************************************************
   ************************************************************
   ************************************************************ */


/* ************************************************************
   ************************************************************
   These functions *DO* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */

int C_midx(int *ii, int *jj, int *nrow);

void C_submat(double *X, int *nrX, int *ncX,
	     int *idx1, int *n1, int *idx2, int *n2, double *ans);

void C_transpose(double *X, int *nrX, int *ncX, double *ans);

void C_symmetrize(double *X, int *nrX);

void C_issym(double *X, int *nrX, int *ncX, double *ans);

void C_matadd(double *X, int *nrX, int *ncX,
	     double *Y, int *nrY, int *ncY, double *ans);

void C_matsubt(double *X, int *nrX, int *ncX,
	      double *Y, int *nrY, int *ncY, double *ans);

void C_schursubt(double *X, int *nrX, int *ncX,
		int *idx1, int *n1, int *idx2, int *n2, double *ans);

/* ************************************************************
   ************************************************************
   These functions *DO NOT* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */

void C_matprod(double *X, int *nrX, int *ncX,
	      double *Y, int *nrY, int *ncY, double *ans);

void C_solve(double *A, int *nrA, double *B, int *ncB);

void C_inverse(double *A, int *nrA);

void C_determinant(double *Ain, int *nrA, double *ans);

void C_traceABsym(double *A, int *nrA, int *ncA,
		 double *B, int *nrB, int *ncB, double *ans);


/* ************************************************************
   ************************************************************
   ************************************************************
   IMPLEMENTATION
   ************************************************************
   ************************************************************
   ************************************************************ */


/* ************************************************************
   ************************************************************
   These functions *DO* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */


/*  midx: */
/*  Maps matrix entry (ii,jj) (running as 1,2,...) from a matrix with */
/*  <nrow> rows into index (running as 0,1,...) */
/* FIXME : DOING from 0...n-1 */

int C_midx(int *ii, int *jj, int *nrow)
{
  //int ans = *ii-1 + *nrow*(*jj-1);
  int ans = *ii + *nrow*(*jj);
  //Rprintf("ii=%d, jj=%d ans=%d\n", *ii, *jj, ans);
  return(ans);
}

SEXP R_midx(SEXP ii, SEXP jj, SEXP nrow)
{

  double *i1 = REAL(ii), *j1=REAL(jj), *nr=REAL(nrow);
  int II = (int) *i1 , JJ = (int) *j1 , NR = (int) *nr ;

  SEXP ans;
  double *ansptr, ANS;

  PROTECT(ans = allocVector(REALSXP, 1));
  ansptr = REAL(ans);

  ANS = (double) C_midx(&II, &JJ, &NR);
  //Rprintf("ANS=%f\n", ANS);

  *ansptr = ANS;

  UNPROTECT(1);
  return(ans);
}


void C_submat(double *X, int *nrX, int *ncX,
	     int *idx1, int *n1, int *idx2, int *n2, double *ans)
{
/*   Rprintf("submat\n"); */
/*   Rprintf("nrX=%d, ncX=%d\n", *nrX, *ncX); */
/*   Rprintf("n1=%d, n2=%d\n", *n1,*n2); */

  //Rprintf("Csubmat\n");
  int ii, jj, kk, idx, ii1,jj1;

  kk = 0;
  for (jj=0; jj < *n2; jj++){
    for (ii=0; ii < *n1; ii++){
      ii1 = (int) idx1[ii];
      jj1 = (int) idx2[jj];
      idx = C_midx(&ii1, &jj1, nrX);
      //      Rprintf("ii=%d, jj=%d idx1=%f idx2=%f, idx=%3d X=%f\n", ii,jj, idx1[ii], idx2[jj], idx, X[idx]);
      ans[kk] = X[idx];
      kk++;
    }
  }
}

SEXP R_submat(SEXP X, SEXP idx1, SEXP idx2)
{

  int *i1, *i2;
  int n1 = length(idx1), n2=length(idx2);
  int    *xdims;
  double *xptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  PROTECT(idx1 = coerceVector(idx1, INTSXP));
  PROTECT(idx2 = coerceVector(idx2, INTSXP));
  i1 = INTEGER(idx1);
  i2 = INTEGER(idx2);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, n1, n2));
  ansptr = REAL(ans);

  C_submat(xptr, &xdims[0], &xdims[1], i1, &n1, i2, &n2, ansptr);

  UNPROTECT(4);
  return(ans);
}


void C_transpose(double *X, int *nrX, int *ncX, double *ans)
{
  //Rprintf("Ctranspose\n");
  int idx1, idx2;
  for (int jj=0; jj< *ncX; jj++){
    for (int ii=0; ii< *ncX; ii++){
      idx1 = C_midx(&ii, &jj, nrX);
      idx2 = C_midx(&jj, &ii, ncX);
      //Rprintf("idx1=%d, idx2=%d\n", idx1, idx2);
      ans[idx2] = X[idx1];
    }
  }
}

SEXP R_transpose(SEXP X)
{

  int    *xdims;
  double *xptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, xdims[1], xdims[0]));
  ansptr = REAL(ans);

  C_transpose(xptr, &xdims[0], &xdims[1], ansptr);

  UNPROTECT(2);
  return(ans);
}

// Symmetrizes X; returns the result X.
void C_symmetrize(double *X, int *nrX)
{
  //Rprintf("Csymmetrize\n");
  int ii;
  double *Xt;
  Xt    = (double *) R_alloc(*nrX**nrX, sizeof(double));
  //Xt    = (double *) Calloc(*nrX**nrX,double);

  C_transpose(X, nrX, nrX, Xt);
  //printmatd(Xt, nrX, nrX);
  for (ii=0; ii < *nrX**nrX; ii++){
    X[ii] = (X[ii] + Xt[ii])/2.0;
  }
  //Free(Xt);
}

SEXP R_symmetrize(SEXP X)
{
  int    ii, *xdims;
  double *xptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, xdims[0], xdims[1]));
  ansptr = REAL(ans);

  for (ii=0; ii<xdims[0]*xdims[0]; ii++){
    ansptr[ii] = xptr[ii];
  }

  C_symmetrize(ansptr, &xdims[0]);

  UNPROTECT(2);
  return(ans);
}


void C_issym(double *X, int *nrX, int *ncX, double *ans)
{
  //Rprintf("Cissym\n");
  int ii;
  double *Xt;
  Xt    = (double *) R_alloc(*nrX**ncX, sizeof(double));
  C_transpose(X, nrX, ncX, Xt);
  *ans = 1;
  for (ii=0; ii<*nrX**ncX; ii++){
    //Rprintf("%f %f\n", X[ii], Xt[ii]);
    if (fabs( X[ii]-Xt[ii] ) > 1e-12){
      //Rprintf(".......\n");
      *ans = 0;
      break;
    }
  }
}

SEXP R_issym(SEXP X)
{
  int    *xdims;
  double *xptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, 1));
  ansptr = REAL(ans);

  C_issym(xptr, &xdims[0], &xdims[1], ansptr);

  UNPROTECT(2);
  return(ans);
}


void C_matadd(double *X, int *nrX, int *ncX,
		double *Y, int *nrY, int *ncY, double *ans)
{
  Rprintf("Cmatadd\n");
  for (int ii=0; ii < *nrX * *ncX; ii++){
    ans[ii] = X[ii]+Y[ii];
  }
}


SEXP Rmatadd(SEXP X, SEXP Y)
{

  int    *xdims, *ydims;
  double *xptr, *yptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  ydims = getDims(Y);
  PROTECT(Y = coerceVector(Y, REALSXP));
  yptr  = REAL(Y);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, xdims[0], xdims[1]));
  ansptr = REAL(ans);

  C_matadd(xptr, &xdims[0], &xdims[1], yptr, &ydims[0], &ydims[1], ansptr);

  UNPROTECT(3);
  return(ans);
}

void C_matsubt(double *X, int *nrX, int *ncX,
		double *Y, int *nrY, int *ncY, double *ans)
{
  Rprintf("Cmatsubt\n");
  for (int ii=0; ii < *nrX * *ncX; ii++){
    ans[ii] = X[ii]-Y[ii];
  }
}

SEXP R_matsubt(SEXP X, SEXP Y)
{
  int    *xdims, *ydims;
  double *xptr, *yptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  ydims = getDims(Y);
  PROTECT(Y = coerceVector(Y, REALSXP));
  yptr  = REAL(Y);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, xdims[0], xdims[1]));
  ansptr = REAL(ans);

  C_matsubt(xptr, &xdims[0], &xdims[1], yptr, &ydims[0], &ydims[1], ansptr);

  UNPROTECT(3);
  return(ans);
}



void C_schursubt(double *X, int *nrX, int *ncX,
	    int *idx1, int *n1, int *idx2, int *n2, double *ans)
{
  //Rprintf("++++++++ in Cschursubt\n");

  double *tmp1, *tmp2;
  double *X22;
  tmp1    = (double *) R_alloc(*n1**n2, sizeof(double));
  tmp2    = (double *) R_alloc(*n1**n2, sizeof(double));
  X22     = (double *) R_alloc(*n2**n2, sizeof(double));
  //double issymX22, issymX22inv, issymX22inv2;

  C_submat(X, nrX, ncX, idx1, n1, idx2, n2, tmp1); // tmp1 = X.12 (n1 x n2)
  //Rprintf("X.12\n"); printmatd(tmp1, n1, n2);

  C_submat(X, nrX, ncX, idx2, n2, idx2, n2, X22);  // X22  = X.22 (n2 x n2)
  //Rprintf("X.22\n"); printmatd(X22, n2, n2);
  //Cissym(X22, n2, n2, &issymX22);

  C_inverse(X22, n2);                              // X22  = inv(X.22) (n2 x n2)
  //Rprintf("inv(X.22)\n"); printmatd(X22, n2, n2);
  //Cissym(X22, n2, n2, &issymX22inv);

  //Csymmetrize(X22, n2);
  //Cissym(X22, n2, n2, &issymX22inv2);
  //Rprintf("issymX22=%f, issymX22inv=%f issymX22inv2=%f\n", issymX22, issymX22inv, issymX22inv2);

  C_matprod(tmp1, n1, n2, X22, n2, n2, tmp2);      // tmp2 = X.12*inv(X.22) (n1 x n2)
  //Rprintf("X.12*inv(X.22)\n"); printmatd(tmp2, n1, n2);

  C_submat(X, nrX, ncX, idx2, n2, idx1, n1, tmp1); // tmp1 = X.21 (n2 x n1)
  //Rprintf("X.21\n"); printmatd(tmp1, n2, n1);

  C_matprod(tmp2, n1, n2, tmp1, n2, n1, ans);      // ans = X.12*inv(X.22)*X.21 (n1 x n1)
  //Rprintf("X.12*inv(X.22)*X.21 \n"); printmatd(ans, n1, n1);

  C_symmetrize(ans, n1);
}

SEXP R_schursubt(SEXP X, SEXP idx1, SEXP idx2)
{

  int *i1, *i2;
  int n1 = length(idx1), n2=length(idx2);

  int    *xdims;
  double *xptr;
  xdims = getDims(X);
  PROTECT(X = coerceVector(X, REALSXP));
  xptr  = REAL(X);

  PROTECT(idx1 = coerceVector(idx1, INTSXP));
  PROTECT(idx2 = coerceVector(idx2, INTSXP));
  i1 = INTEGER(idx1);
  i2 = INTEGER(idx2);

  double *ansptr;
  SEXP ans;

  PROTECT(ans = allocMatrix(REALSXP, n1, n1));
  ansptr = REAL(ans);

  C_schursubt(xptr, &xdims[0], &xdims[1], i1, &n1, i2, &n2, ansptr);

  UNPROTECT(4);
  return(ans);
}



/* ************************************************************
   ************************************************************
   These functions *DO NOT* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */

/* Cmatprod: */
/* Computes X * Y; returns answer in ans */
void C_matprod(double *X, int *nrX, int *ncX,
		double *Y, int *nrY, int *ncY, double *ans)
{
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, nrX, ncY, ncX, &one,
		  X, nrX, Y, nrY, &zero, ans, nrX FCONE FCONE);
}

/* Csolve: */
/* Solves the system Ax=B; returns the answer in B */
void C_solve(double *A, int *nrA, double *B, int *ncB)
{
  //Rprintf("Csolve\n");
  int info;
  int *ipiv  = (int *) R_alloc((int) *nrA, sizeof(int));
  //int *ipiv  = (int *) malloc (*nrA * sizeof(int));
  F77_CALL(dgesv)(nrA, ncB, A, nrA, ipiv, B, nrA, &info);
  //free(ipiv);
}

/* Cinverse: */
/* Finds the inverse of A; returns the answer in A */
void C_inverse(double *A, int *nrA)
{
  //Rprintf("Cinverse\n");
  int ii, jj, info, *ipiv, n=*nrA, n2=*nrA**nrA;
  double *B;
  B     = (double *) R_alloc(n2, sizeof(double)); //Rprintf("B ok...\n");
  ipiv  = (int *)    R_alloc(n, sizeof(int));     //Rprintf("ipiv ok...\n");

  for (ii=0; ii<*nrA; ii++){
    for (jj=0; jj<*nrA; jj++){
      if (ii==jj){
	B[ii+*nrA*jj] = 1;
      } else {
	B[ii+*nrA*jj] = 0;
      }
    }
  }

  F77_CALL(dgesv)(nrA, nrA, A, nrA, ipiv, B, nrA, &info);
  //Memcpy(A, B, (size_t) n2);
  for (ii=0; ii < n2; ii++)
    A[ii] = B[ii];
  //Rprintf("Cinverse-done\n");
}


void C_determinant(double *Ain, int *nrA, double *ans)
{
  int i, n, info, *jpvt, sign, useLog=1;
  double modulus = 0.0; /* -Wall */
  int nrA2;
  double *Awork;

  nrA2 = *nrA * *nrA;

  Awork = (double *) R_alloc(*nrA * *nrA, sizeof(double));
  Memcpy(Awork, Ain, (size_t) (nrA2));

  n = *nrA;
  //jpvt = (int *) R_alloc(*nrA, sizeof(int));

  jpvt = (int *) malloc(*nrA *  sizeof(int));
  F77_CALL(dgetrf)(&n, &n, Awork, &n, jpvt, &info);

  //printmatd(Awork, &n, &n);

  sign = 1;
  if (info < 0)
    Rprintf("error code %d from Lapack routine '%s'", info, "dgetrf");
  else if (info > 0) { /* Singular matrix:  U[i,i] (i := info) is 0 */
    /*warning("Lapack dgetrf(): singular matrix: U[%d,%d]=0", info,info);*/
    modulus = (useLog ? R_NegInf : 0.);
  }
  else {
    for (i = 0; i < n; i++) if (jpvt[i] != (i + 1))
      sign = -sign;
    if (useLog) {
      modulus = 0.0;
      for (i = 0; i < n; i++) {
	double dii = Awork[i*(n + 1)]; /* ith diagonal element */
	modulus += log(dii < 0 ? -dii : dii);
	if (dii < 0) sign = -sign;
      }
    } else {
      modulus = 1.0;
      for (i = 0; i < n; i++)
	modulus *= Awork[i*(n + 1)];
      if (modulus < 0) {
	modulus = -modulus;
	sign = -sign;
      }
    }
  }

  //Rprintf("%i %f %f \n", sign, modulus, ((float) sign) * exp(modulus));
  *ans = sign * exp(modulus);
  free(jpvt);
}


// Calculates tr(AB) where A and B *must* be symmetrical
void C_traceABsym(double *A, int *nrA, int *ncA,
		 double *B, int *nrB, int *ncB, double *ans)
{
  int ii;
  double x=0;
  
  for (ii=0; ii<*nrA**nrA; ii++){
    x=x+A[ii]*B[ii];
  }

  //Rprintf("tr: %f\n", x);
  *ans=x;
}
