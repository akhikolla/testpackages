#include <Rdefines.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <R_ext/Lapack.h> 
//#include <R_ext/Applic.h>       /* for dgemm */
//#include <R_ext/RS.h>		/* for F77_... */

#include "_utils_mat.h"
#include "_utils_print.h"


double eval_logL(double *S, double *K, int *nobs, int *nvar, int *details);

double eval_maxparmdiff(double *currKinv, double *prevKinv, int *nvar);

void update_K(double *S, double *K, int *nvar, int *nobs,
	     int *gset, int *ng, int *cset, int *nc,
	     double *Swork, double *Dwork, int *details);

void Cggmfit(double *S,     int *nobs,    double *K, 
	     int *nvar,      int *ngen, 
	     int *glen,      int *glist, 
	     int *clen,      int *clist, 
	     double *logL,   double *eps,  int *iter, 
	     int *converged, int *details)
{
  // Rprintf("WORKING VERSION OF CGGMFIT\n");

  int ii, kk, currg, ng, nc, itcount=0, firstg, lastg, firstc, lastc;
  int *gset, *cset; 
  int nvar2 = *nvar * *nvar;
  double *Dwork, *Swork, prevlogL;
  double *prevKinv, *currKinv;
  double maxparmdiff;

  /*   Rprintf("OXFORD: S (start):\n"); printmatd(S, nvar, nvar); */
  /*   Rprintf("OXFORD: K (start):\n"); printmatd(K, nvar, nvar); */

  // Holding generator and complement
  gset  = (int *) R_alloc(*nvar, sizeof(int));
  cset  = (int *) R_alloc(*nvar, sizeof(int));

  // Holding previous and current K (used for maxparmdiff)
  prevKinv = (double*) R_alloc (nvar2, sizeof(double));
  currKinv = (double*) R_alloc (nvar2, sizeof(double));

  // Temporary matrices
  Dwork = (double*) R_alloc (nvar2, sizeof(double));
  Swork = (double*) R_alloc (nvar2, sizeof(double));

  // Initial log-likelihood
  prevlogL = eval_logL(S, K, nobs, nvar, details);
  if (*details>=1)
    Rprintf(". Initial logL: %14.6f \n", prevlogL);

  if (*ngen==1) // It is the saturated model
    {
      for (ii=0; ii<*nvar * *nvar; ii++){
	K[ii] = S[ii];
      }
      C_inverse(K, nvar);
      *logL = eval_logL(S, K, nobs, nvar, details);
    }
  else
    {
      Memcpy(prevKinv, K, (size_t) nvar2);
      C_inverse(prevKinv, nvar);

      while(1){
	itcount++;
	firstg = 0;
	firstc = 0;
	lastg  = 0;
	lastc  = 0;
	
	for (currg = 0; currg < *ngen; currg++)
	  {
	    ng    = glen[currg]; // No. of elements in generator
	    nc    = clen[currg]; // No. of elements in complement
	    
	    // Find 'positions' of current generator (in glist) and complements (in clist) 
	    firstg = lastg;
	    firstc = lastc;
	    lastg  = lastg+glen[currg];
	    lastc  = lastc+clen[currg];
	    
	    // Copy indices of current generator and current complement to gset and cset
	    kk = 0;
	    for (ii=firstg; ii<lastg; ii++){
	      gset[kk++] =  glist[ii];
	    }  
	    kk = 0;
	    for (ii=firstc; ii<lastc; ii++){
	      cset[kk++] =  clist[ii];
	    } 
	    
	    // Update K
	    update_K(S,K,nvar,nobs, gset,&ng,cset,&nc, Swork, Dwork, details);
	  } /* for */

	// Finding max parameter difference
	Memcpy(currKinv, K, (size_t) nvar2);
	C_inverse(currKinv, nvar);
	maxparmdiff = eval_maxparmdiff(currKinv, prevKinv, nvar);
	Memcpy(prevKinv, currKinv, (size_t) nvar2);
	
	// Find log-likelihood
	*logL = eval_logL(S, K, nobs, nvar, details );
	
	//Rprintf("det %f logdet %f tr %f logL %f n %i \n", det, log(det), trAB, *logL, *n);
	if (*details>=1)
	  Rprintf(". Iteration: %3i logL: %14.6f diff logL: %20.13f maxparmdiff: %18.12f\n", 
		  itcount, *logL, *logL-prevlogL, maxparmdiff);
	
	if (itcount>=1){
	    if (maxparmdiff < *eps){
	    *converged = 1;
	    break;
	  } else {
	    if (itcount==*iter){
	      *converged = 0;
	      break;
	    }
	  }      
	}
	prevlogL = *logL;
      }
      if (*details>=1)
	Rprintf(". Final: Iterations: %i logL: %f diff logL: %16.13f  maxparmdiff: %18.12f\n", 
		itcount, *logL, *logL-prevlogL, maxparmdiff);
    }
  *iter = itcount;
}



double eval_logL(double *S, double *K, int *nobs, int *nvar, int *details)
{

  
  double det, trAB, logL;
  
  C_determinant(K, nvar, &det);
  C_traceABsym(K, nvar, nvar, S, nvar, nvar, &trAB);
  double diff = (double) log(det) - trAB;

  double nobs22 = (double) *nobs / 2;
  
  //Rprintf("nobs=%d nvar=%d\n", *nobs, *nvar);
  logL =  - (*nobs**nvar*log(2*M_PI)/2) + nobs22 * diff;
  //logL =  nobs22 * diff;
  if (*details>=3){
    Rprintf(" ...in eval_logL:\n");
    Rprintf(" ...in eval_logL: det=%16.12f, logdet=%16.12f, tr=%16.12f diff=%16.12f, nobs=%d, logL=%f\n", 
	    det, log(det), trAB, diff, *nobs, logL);
    /*   Rprintf("in eval_logL: S :\n"); printmatd(S, nvar, nvar); */
    /*   Rprintf("in eval_logL: K :\n"); printmatd(K, nvar, nvar); */
  }
  return(logL);
}

void update_K(double *S, double *K, int *nvar, int *nobs,
	     int *gset, int *ng, int *cset, int *nc,
	     double *Swork, double *Dwork, int *details)
{
  
  //Rprintf("---- in update_K ----\n");
  int ii, jj, kk;

  // S(gg) -> Swork
  C_submat(S, nvar, nvar, gset, ng, gset, ng, Swork);
  //Rprintf("(Sgg):\n"); printmatd(Swork, ng, ng);
  
  // inv(S(gg)) -> Swork
  C_inverse(Swork, ng);
  //Rprintf("inv(S(gg)):\n"); printmatd(Swork, ng, ng);
	
  // K(gc)inv(K(cc))K(cc) -> Dwork
  C_schursubt(K, nvar, nvar, gset, ng,  cset, nc, Dwork);
  //Rprintf("schur:\n"); printmatd(Dwork, ng, ng); 

  // Update: K(gg) <- inv(S(gg)) + K(gc)inv(K(cc))K(cc)
  kk = 0;
  for (jj=0; jj<*ng; jj++){
    for (ii=0; ii<*ng; ii++){
      //Rprintf(" %i %i %i \n", gset[ii], gset[jj], kk);
      K[(int)(gset[ii] + *nvar * gset[jj])] = Dwork[kk] + Swork[kk];
      kk++;
    }
  }
  
  if (*details>=2){
    double innerlogL = eval_logL(S, K, nobs, nvar, details); 
    Rprintf(".. updating generator :");
    printveci(gset, ng); Rprintf(" //"); printveci(cset, nc); 
    Rprintf("logL (after update)=%16.12f\n", innerlogL);
  }
}

double eval_maxparmdiff(double *currKinv, double *prevKinv, int *nvar)
{
  double maxparmdiff, absdiff;
  int ii, jj;
  int ii_idx, jj_idx, ij_idx;
  maxparmdiff = 0;
  for (jj=0; jj< *nvar; jj++){
    for (ii=0; ii< *nvar; ii++){
      ij_idx = C_midx(&ii, &jj, nvar);
      ii_idx = C_midx(&ii, &ii, nvar);
      jj_idx = C_midx(&jj, &jj, nvar);
      //Rprintf("ij_idx=%3d ii_idx=%3d, jj_idx=%3d\n", ij_idx, ii_idx, jj_idx);
      double num = fabs(currKinv[ij_idx]-prevKinv[ij_idx]);
      double den = sqrt(currKinv[ii_idx]*currKinv[jj_idx]+currKinv[ij_idx]*currKinv[ij_idx]);
      absdiff = num/den;
      //Rprintf("num=%18.8f den=%18.8f, absdiff=%18.8f ij=%5f ii=%f jj=%f\n", 
      //    num, den, absdiff, currKinv[ij_idx], currKinv[ii_idx], currKinv[jj_idx]);
      if( absdiff > maxparmdiff ){
	maxparmdiff = absdiff;
      }	    
    }
  }
  return(maxparmdiff);
}


  //Rprintf("update_K: K (after update):\n"); printmatd(K, nvar, nvar);
/*   double issymK=0, issym_schur=0; */
/*   Cissym(K, nvar, nvar, &issymK);  */
/*   Cissym(Dwork, ng, ng, &issym_schur); */
  // Rprintf("---- update_K (exit) issymK=%f issym_schur=%f\n", issymK, issym_schur);


  /*   Find sizes of temporary vectors  */
  /*    We don't use that anymore... */
  /*   for (ii=0; ii<*ngen; ii++){ */
  /*     xx = glen[ii] * glen[ii];    //Rprintf(" xx %i\n",xx); */
  /*     if (xx > maxn_gg) */
  /*       maxn_gg = xx; */
  /*     xx = clen[ii] * clen[ii];    //Rprintf(" xx %i\n",xx); */
  /*     if (xx > maxn_cc) */
  /*       maxn_cc = xx; */
  /*     xx = glen[ii] * clen[ii];    //Rprintf(" xx %i\n",xx); */
  /*     if (xx > maxn_gc)  */
  /*       maxn_gc = xx;      */
  /*   }  */
  /*   Dwork = (double*) R_alloc (maxn_gg, sizeof(double)); */
  /*   Swork = (double*) R_alloc (maxn_gg, sizeof(double)); */








































