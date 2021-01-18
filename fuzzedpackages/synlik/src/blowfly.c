/* This file contains a model definition file for a stochastic Nisbet and Gurney
blowfly model.
The idea is to efficiently solve the  model for multiple
replicates, conditional on an unscaled noise vector.
*/
  #include <R.h>
  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include "summ_stats.h"
  
void blowC(double *n,double *theta,double *e,double *e1,int *burn_in,int *n_t, int *n_reps) {
/* Simulates `n_reps', length `n_t' replicates of blowfly model, discarding an initial 
   sequence of length `burn_in'.
   e (and e1, if used) must be length (burn_in+n_t)*n_reps: these are the noise sequences - typically 
   gamma random variables with mean 1, but variance set externally. 

   n is n_t by n_reps
   This should be much faster than looping in R, even if all reps parallel.
   In this version the noise is an AR1 proccess --- the raw error stream is
   E_{t+1} = alpha E_t + e_t * (1 - alpha), where 0 <= alpha < 1. E_t is scaled 
   by sig before use. 
  
*/
  int i,j,lag;
  double delta,P,N_0,sig,tau,*N,*N1,*N0,*Nlag,sig1;
  delta = theta[0];
  P = theta[1];
  N_0 = theta[2];
  sig = theta[3];
  tau = theta[4]; /* note that this is discrete */
  sig1 = theta[5]; /* the second noise parameter --- used externally */
  lag = (int) floor(tau);if (tau-lag>.5) lag++;
  /* Allocate single replicate storage */
  N = (double *)calloc((size_t)lag + *burn_in + *n_t,sizeof(double));
  for (i=0;i<lag;i++) N[i] = 180.0;
  /* following iterates the blowfly model... */
  for (j=0;j<*n_reps;j++) {
    N1 = N + lag;N0 = N1 - 1;Nlag=N;
    for (i=lag;i< *burn_in+lag;i++,N0++,N1++,Nlag++,e++,e1++) {
      /* This is old discretization, less stable, and less natural than one used now...
         *N1 = *N0 * exp(P * *Nlag * exp(- *Nlag/ N_0)/ *N0 - delta + err * sig); */
      
      *N1 = P * *Nlag * exp(- *Nlag/ N_0) * *e + *N0 * exp( - delta * *e1 );
  
    }
    for (i=0;i<*n_t;i++,n++,N0++,N1++,Nlag++,e++,e1++) {
      /* *n =  *N1 = *N0 * exp(P* *Nlag*exp(- *Nlag/N_0)/ *N0 - delta + err * sig); */
      *n =  *N1 =  P * *Nlag * exp(- *Nlag/ N_0) * *e + *N0 * exp( - delta * *e1);
    }
  }
  free(N);
} /* end of ng_bf */


