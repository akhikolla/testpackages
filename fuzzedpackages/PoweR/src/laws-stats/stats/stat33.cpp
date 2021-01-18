// Title: Statistique de test de Gel-Miao-Gastwirth
// Ref. (book or article): Gel, Y.R., Miao, W. and Gastwirth, J.L. (2007), Robust directed tests of normality against heavy-tailed alternatives, 
//						   Computational Statistics & Data Analysis, Vol. 51, pp. 2734-2746.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat33(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
   alter[0] = 3;

   int i, j=0, n=xlen[0];
   if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
     const char *nom = "$R_{sJ}$";
// Here, INDICATE the number of parameters of your statistic
     nbparamstat[0] = 0;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).

     }
// The following 7 lines should NOT be modified
     const char *space = " ";
     while (nom[j] != '\0') {
       name[j][0] = nom[j];
       j++;
     }
     for (i=j;i<50;i++) name[i][0] = space[0];
     return;
   }

    if (n>3) {
// Computation of the value of the test statistic
      void R_rsort (double* x, int n);
      //      double pnorm(double q, double mean, double sd, int lower_tail, int log_p);      
      double *xtmp;
      xtmp = new double[n];
      double statRsJ, varX=0.0, meanX=0.0, Jn=0.0, M, tmp, sdX, pi;
	  pi = 4.0*atan(1.0); 		// or use pi = M_PI, where M_PI is defined in math.h

	  // calculate sample mean
      for (i=0;i<=(n-1);i++) {
	    meanX = meanX + x[i];
      }
	  meanX =meanX/(double)n;
      
	  // calculate sample var and standard deviation
	  // for (i=0;i<=(n-1);i++) {
	    // varX = varX + R_pow(x[i],2.0);
      // }
	  // varX = ((double)n)*(varX/(double)n - R_pow(meanX,2.0))/(double)(n-1); 
      // sdX = sqrt(varX);
	  for (i=0;i<n;i++) {
	    varX = varX + R_pow(x[i]-meanX,2.0);
	  }
	  varX = varX/(double)n;	
	  sdX = sqrt(varX);
	  
	  // calculate sample median
      for (i=0;i<=(n-1);i++) {
	    xtmp[i] = x[i];
      }
	  
      R_rsort(xtmp,n); // We sort the data
     
	  if ((n%2) == 0) {
	    M = (xtmp[n/2]+xtmp[n/2-1])/2.0; 
	  } else {
	    M = xtmp[n/2]; // sample median
      }
	  
	  // calculate statRsJ
	  for (i=0;i<=(n-1);i++) {
	    Jn = Jn + fabs(x[i]-M);
      }
	  Jn = sqrt(pi/2.0)*Jn/((double)n);
      
	  statRsJ = sdX/Jn;
      
      statistic[0] = statRsJ; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue33.cpp"
}
    
// We take the decision to reject or not to reject the null hypothesis H0
    // for (i=0;i<=(nblevel[0]-1);i++) {
     
    // if (usecrit[i] == 1) { // We use the provided critical values
	// if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	// } else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	// } else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
      // } else {
	    // if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      // }
	  	  
    // }
	
// We take the decision to reject or not to reject the null hypothesis H0
     for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }

 // If applicable, we free the unused array of pointers
      delete[] xtmp;

    }
    
// We return
    return;
   
    
  }
  
}
