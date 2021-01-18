// Title: The Gulati statistic for the Laplace distribution
// Ref. (book or article): Gulati, Sneh (2011), Goodness of fit test for the Rayleigh and the Laplace distributions,
//						   International Journal of Applied Mathematics & Statistics, Vol. 24, No. SI-11A, pp. 74-85. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat59(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$Z$";
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
      double *Y;
      Y = new double[n];
      double *W;
      W = new double[n];
      double *U;
      U = new double[n-1];
      double statZ, muhat, bhat, tmp = 0.0, tn = 0.0, ubar, tmp2 = 0.0, tmp3 = 0.0, Z1square, Z2square, tmp4 = 0.0;
      
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
	
      // generate vector Y
      for (i = 0; i < n; i++) {
	Y[i] = fabs(x[i] - muhat) / bhat;
      }
      R_rsort(Y, n); 			// we sort the data 
      
      // generate vector W
      W[0] = (double)n *(Y[0] - 0.0);
      tn = W[0];
      for (i = 1; i < n; i++) {
	W[i] = ((double)n - (double)i) * (Y[i] - Y[i - 1]);
	tn = tn + W[i];
      }
      
      for (i = 0; i < (n - 1); i++) {
	tmp2 = tmp2 + W[i];
	U[i] = tmp2 / tn;
	tmp3 = tmp3 + U[i];
      }
      ubar = tmp3 / (double)(n - 1);
	
      // calculate statZ
      Z1square = (double)(12 * (n - 1)) * R_pow(ubar - 0.5, 2.0);
      
      for (i = 0; i < (n-1); i++) {
	tmp4 = tmp4 + (double)(i + 1) * U[i];
      }
      Z2square = ((double)(5 * (n - 1)) / (double)((n + 1) * (n - 2))) * R_pow((double)n - 2.0 + 6.0 * (double)n * ubar - 12.0 * tmp4 / (double)(n - 1), 2.0);
      
      statZ = Z1square + Z2square;
      
      statistic[0] = statZ; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue59.cpp"
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
    delete[] Y;
    delete[] W;
    delete[] U;

}

// We return
    return;
   
        
  }
  
  
  
}
