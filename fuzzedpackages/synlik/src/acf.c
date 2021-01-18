#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "summ_stats.h"

/* Routine for vectorized ACF calculation */

void slacf(double *acf,double *x,int *n,int *n_reps,int *max_lag,double *NAcode,
           int *correlation)
/* `x' is a n by n_reps matrix, stored column wise.
       Each column of x is one series. 
      
   `max_lag' gives the maximum lag to calculate the correlation for.
   
   `NAcode' is the number standing for NA. All NA series will cause failure.

   `correlation' should be 0 for auto-covariance, 1 for auto-correlation

    `acf' contains the acf coefficients, stored columnwise.
*/
{ double *xm,*p,xx,*p0,*p1,*p2;
  int i,j,k,lag;
  xm = (double *)calloc((size_t)*n_reps,sizeof(double)); /* storage for col means */

  /* get column (series) means */
  for (p=x,j=0;j<*n_reps;j++) {
    for (k=0,xx=0.0,i=0;i<*n;i++,p++) if (*p != *NAcode) { xx += *p;k++;}
    xm[j] = xx/k;  
  }

  /* centre columns */
  for (p=x,j=0;j<*n_reps;j++) {
    for (xx=xm[j],i=0;i<*n;i++,p++) if (*p != *NAcode) *p -= xx;
  }
   
  /* get covariances */
  
  for (p=x,i=0;i<*n_reps;i++,p+=*n) /* loop through series */ 
    for (lag=0;lag<=*max_lag;lag++) 
      { for (k=0,xx=0.0,p0=p,p1=p + *n - lag,p2=p+lag;p0<p1;p0++,p2++) 
          if (*p0!=*NAcode&&*p2!=*NAcode) { k++;xx += *p0 * *p2;}
	if (k>0) acf[lag + (*max_lag+1) * i] = xx/k;
	  else  acf[lag + (*max_lag+1) * i] = *NAcode;
      }
 
  /* convert to correlations? */
  if (*correlation) 
  for (p=acf,i=0;i<*n_reps;i++)
    for (xx = *p,p1=p + *max_lag;p<=p1;p++) *p/=xx;

  free(xm);   
}


