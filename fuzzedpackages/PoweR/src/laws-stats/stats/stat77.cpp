// Title: The Swartz statistic for uniformity
// Ref. (book or article): Swartz, T. (1992), Goodness-of-fit tests using Kullback-Leibler information,
//                         Communications in Statistics. Theory and Methods, 21, 711-729.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat77(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$A^{*}(n)$";
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
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *U;
	double *Z1;
    double *Z2;
	U = new double[n];
	Z1 = new double[n];
	Z2 = new double[n];
	double statAn, sumAn=0.0;
	

	// generate vector U
	for (i=0;i<n;i++) {
	  U[i] = punif(x[i],0.0,1.0,1,0);
	}
    R_rsort(U,n); // We sort the data
	
	    
	// calculate statAn
	for (i=0; i<n; i++) {
	  if (i < (n-1)) { 
	    Z1[i] = U[i+1];
	  } else Z1[i] = 2.0 - U[n-1];
	  if (i > 0) {
	    Z2[i] = U[i-1];
	  } else Z2[i] = -U[0];
	  
	  sumAn = sumAn + R_pow((Z1[i]-Z2[i])/2.0 - 1.0/(double)n,2.0);
    
	}

	statAn = (double)n*sumAn/2.0;
	
    statistic[0] = statAn; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue77.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
        }
    }
    
// If applicable, we free the unused array of pointers
    delete[] U;
	delete[] Z1;
	delete[] Z2;

}

// We return
    return;
   
        
  }
  
  
  
}
