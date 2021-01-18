#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "summ_stats.h"

/* Code for polynomial auto-regression*/

void slnlar(double *beta, double *x,int *n,int *n_reps,int *n_terms,
            int *lag,int *power,double *NAcode)
/* `x' is a n by n_reps matrix, stored column wise.
       Each column of x is one series. 
      
   `n_terms' gives the number of terms on the rhs of the autoregression.

   `lag[i]' gives the lag of the ith term on the rhs.

   `power[i]' gives the power to which the ith rhs term should be raised.
   
   `NAcode' is the number standing for NA. All NA series will cause failure.

   `beta' contains the ar coefficients, stored columnwise.
*/
{ int i,j,k,nx,max_lag,ny,*na,ip,*pi,*pi1,one=1,*pivot,ok;
  double *X,*y,*xp,*Xp,xx,*p,*p1,*p2,*b,*tau,obs1;

  for (max_lag=0,i=0;i<*n_terms;i++) if (lag[i]>max_lag) max_lag=lag[i];
  ny = *n - max_lag; /* maximum response vector length */

  Xp = X = (double *)calloc((size_t) *n_terms * ny,sizeof(double)); /* model matrix */
  na = (int *)calloc((size_t) ny, sizeof(int)); 
  tau = (double *)calloc((size_t)*n_terms,sizeof(double));
  pivot = (int *)calloc((size_t)*n_terms,sizeof(int));
  b = (double *)calloc((size_t)*n_terms,sizeof(double));


  for (y=x+max_lag,k=0;k<*n_reps;k++,y += *n) { /* loop through the series */
    /* centre the whole series (y is just the response vector) */
    for (i=0,xx=0.0,p=y-max_lag,p1=p + *n;p<p1;p++) if (*p != *NAcode) { i++;xx += *p;}
    xx/=i; /* series mean */

    ok = 0;obs1 = *NAcode;
    for (p=y-max_lag,p1=p + *n;p<p1;p++) if (*p!=*NAcode) { 
      if (!ok&&p<p1-max_lag) { if (obs1 == *NAcode) obs1=*p; else if (*p!=obs1) ok=1;} /* check that data not all same */
      *p += -xx; /* subtracting series mean */
    } 

    if (ok) { /* data not all the same */
      /* clear the na field */
      for (pi = na,pi1 = na + ny;pi < pi1;pi++) *pi = 0;

      /* build the model matrix, noting NA rows */
      for (Xp=X,i=0;i<*n_terms;i++) { /* work through the terms */
        xp = y - lag[i]; /* the current predictor */
        for (j=0;j<ny;j++,Xp++,xp++) {
          *Xp = 1.0;
          if (*xp == *NAcode) na[j] = 1; else {
            for (ip=0;ip<power[i];ip++) *Xp *= *xp;
          }
        }
      }

      /* add response NAs to na field */
      for (p=y,pi = na,pi1 = na + ny;pi < pi1;pi++,p++) 
	if (*p == *NAcode) *pi = 1; /* set na[i] to 1 is y[i] is NAcode */
      
      /* drop the NA rows from the model matrix, and response data*/
      for (nx=0,pi=na,p=y,p1=y,p2=y+ny;p<p2;p++,pi++) 
        if (!*pi) {*p1=*p;p1++;nx++;} /* drop NAs from y */
      
      for (Xp=X,p=X,i=0;i<*n_terms;i++)
        for (pi=na,j=0;j<ny;j++,pi++,p++) if (!*pi) {*Xp=*p;Xp++;} /* drop NAs from X */
    
      /* QR decompose the model matrix */
   
      mgcv_qr(X,&nx,n_terms,pivot,tau);
      /* solve R b = Q'y for b */
      mgcv_qrqy(y,X,tau,&nx,&one,n_terms,&one,&one); /* y <- Q'y */
      mgcv_backsolve(X,&nx,n_terms,y,b,&one); /* b <- R^{-1} Q'y */

      /* store b in kth column of beta */
      for (i=0;i<*n_terms;i++) beta[k * *n_terms + pivot[i]] = b[i]; 
    } else { /* data had no variability */
      for (i=0;i<*n_terms;i++) beta[k * *n_terms + i] = 0.0;

    }

  }
  free(X);free(tau);free(na);free(pivot);free(b);
} /* end of slnlar */



int d_cmp(const void *p1, const void *p2)
{ /* double comparison function for qsort */
   if (*(double*)p1 < *(double*)p2) return(-1);
   else
   if (*(double*)p1 > *(double*)p2) return(1);
   else return(0);
}


void order_reg(double *beta, double *x,double *z,int *n,int *n_reps,int *np,int *diff) {
/* 1. x is an n by n_reps matrix, where each column is a replicate series.
      z is an n vector. x and z contain no NA's
   2. This routine first differences z and the columns of x `diff' times. 
   3. z and each column of x are then centred and sorted. 
   4. Each column of x is then regressed on z using a polynomial regression of 
      order np. 
   5. A matrix of regression coefficients is returned in np by n_reps matrix beta 
   The idea is that the regression coefficients summarize the shape of the cdf of 
   the differenced series.
*/
  int i,j,k,d,one=1,*pivot;
  double *p,*p1,*p2,*p3,xx,*X,*tau;
  k = *n;
  for (d=0;d<*diff;d++) { /* differencing loop */
    for (p=z,p1=z+1,p2=z+k;p1<p2;p++,p1++) *p = *p1 - *p; /* i.e. z[i+1]-z[i] */
    for (p=x,i=0;i<*n_reps;i++) {
      for (p1=x+i*k,p2=p1+1,p3=p1+k;p2<p3;p1++,p2++,p++) *p = *p2 - *p1;
    }
    k--;
  }
  /* end of differencing. k is now the number of differences */

  for (xx=0.0,p=z,p1=z+k;p<p1;p++) xx += *p;
  xx/=k;for (p=z,p1=z+k;p<p1;p++) *p -= xx; /* z centred */
  for (i=0;i<*n_reps;i++) { /* centre the columns of x */
     for (xx=0.0,p=x+i*k,p1=p+k;p<p1;p++) xx += *p;
     xx/=k;for (p=x+i*k,p1=p+k;p<p1;p++) *p -= xx;
  }

  /* centring complete, now sort... */
  qsort(z,(size_t)k,sizeof(double),d_cmp);
  for (p=x,i=0;i<*n_reps;i++,p+=k) qsort(p,(size_t)k,sizeof(double),d_cmp);

  /* DEBUG CODE... */
  //printf("\n");
  //for (i=0;i<k;i++) printf("%g  ",x[k+i]);
    //  printf("\n");


  /* Now create the model matrix, using contents of `z' */
  if (*np<1) *np = 1;
  X = (double *)calloc((size_t)*np * k,sizeof(double));
  for (p=X,p1=z,p2=z+k;p1<p2;p++,p1++) *p = *p1;
  for (p=X,p1=X+k,i=1;i<*np;i++) 
    for (p2=z,p3=z+k;p2<p3;p++,p1++,p2++) *p1 = *p * *p2;

  /* QR decompose the model matrix */
  tau = (double *)calloc((size_t)*np,sizeof(double));
  pivot = (int *)calloc((size_t)*np,sizeof(int));
  mgcv_qr(X,&k,np,pivot,tau);
  
  /* solve R b = Q'x for b */
  mgcv_qrqy(x,X,tau,&k,n_reps,np,&one,&one); /* y <- Q'y */

  /* need to discard the last k-np rows of x... */
  for (p1=x,p=x,i=0;i<*n_reps;i++,p+= k - *np) 
    for (j=0;j<*np;j++,p1++,p++) *p1 = *p;
 
  mgcv_backsolve(X,&k,np,x,beta,n_reps); /* b <- R^{-1} Q'y */

  /* unpivot the columns of beta */
    for (p=beta,j=0;j<*n_reps;j++,p+=*np) {
    for (i=0;i<*np;i++) tau[i] = p[pivot[i]];
    for (i=0;i<*np;i++) p[i] = tau[i];
  } 
  free(X);free(tau);free(pivot);
} /* end of order_reg */
