// Title: The Gel statistic for the Laplace distribution
// Ref. (book or article): Gel, Yulia R. (2010), Test of fit for a Laplace distribution against heavier tailed alternatives,
//						   Computational Statistics & Data Analysis, Vol. 54, No. 4, 958--965. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat60(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$K$";
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
      for (i = j; i < 50; i++) name[i][0] = space[0];
      return;
    }

    if (n > 3) {
// Computation of the value of the test statistic
      void R_rsort (double* x, int n);
      double statK, muhat, bhat, zeta, xbar, v1sqrt, v2, tmp=0.0, tmp2=0.0, tmp3=0.0, tmp4=0.0, C1, C2;

      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^
      R_rsort(x, n); 		// we sort the data
      if (n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[n / 2 - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[n / 2];
      }
      
      // calculate b^
      for (i = 0; i < n; i++) {
	tmp = tmp + fabs(x[i] - muhat);
      }
      bhat = tmp / (double)n;
	
      // calculate zeta = ML estimate of population standard deviation
      zeta = sqrt(2.0) * bhat;
      
      // calculate xbar
      for (i = 0; i < n; i++) {
	tmp2 = tmp2 + x[i];
      }
      xbar = tmp2 / (double)n;
      
      // calculate v1sqrt and v2
      for (i = 0; i < n; i++) {
	tmp3 = tmp3 + R_pow(x[i] - xbar, 3.0);
	tmp4 = tmp4 + R_pow(x[i] - xbar, 4.0);
      }
      
      v1sqrt = (tmp3 / ((double)n) * R_pow(zeta, 3.0));
      v2 = (tmp4 / ((double)n) * R_pow(zeta, 4.0));
      
      // calculate statK
      C1 = 60.0;
      C2 = 1200.0;
      
      statK = (double)n * R_pow(v1sqrt, 2.0) / C1 + (double)n * R_pow(v2 - 6.0, 2.0) / C2;	
      
      statistic[0] = statK; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue60.cpp"
      }

// We take the decision to reject or not to reject the null hypothesis H0
      for (i = 0; i < nblevel[0]; i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
	} else {
	  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
        }
      }
    
// If applicable, we free the unused array of pointers

}

// We return
    return;
    
  }
  
}
