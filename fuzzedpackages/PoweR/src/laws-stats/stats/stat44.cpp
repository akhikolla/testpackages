// Title: The Watson statistic for the Laplace distribution
// Ref. (book or article): Puig, P. and Stephens, M. A. (2000). Tests of fit for the Laplace distribution, with applications. Technometrics 42, 417-424.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat44(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$U^2$";
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
      double plaplace(double y);
      double *Wa;
      Wa = new double[n];
      double statW, tmp = 0.0, bhat, muhat, sumW = 0.0, meanW = 0.0;
      
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
	
      // generate vector Wa
      for (i = 0; i < n; i++) {
	Wa[i] = plaplace((x[i] - muhat) / bhat);
      }
      
      // calculate statW
      for (i = 1; i <= n; i++) {
	sumW = sumW + R_pow((double)(2 * i - 1) / (double)(2 * n) - Wa[i-1], 2.0);
	meanW = meanW + Wa[i - 1];
      }
      meanW = meanW / double(n);
      
      statW = 1.0 / (double)(12 * n) + sumW - (double)n * R_pow(meanW - 0.5, 2.0);
      
      statistic[0] = statW; // Here is the test statistic value
      
      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue44.cpp"
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
      delete[] Wa;
      
    }
    
    // We return
    return;
    
        
  }
  
  
  // In stat42.cpp, we already defined this function so no need to include it here.
  // The cumulative Laplace distribution function with \mu = 0 and \theta = 1
  // fabs = returns the absolute value of x (a negative value becomes positive, positive value is unchanged). 
  // double plaplace(double x) {
    // double temp = 0.5 * exp(-fabs(x));
    // return (x <= 0.0) ? temp : 1.0 - temp;
  // }
  
  
}
