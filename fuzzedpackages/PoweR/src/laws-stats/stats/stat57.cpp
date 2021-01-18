// Title: The Langholz-Kronmal statistic for the Laplace distribution
// Ref. (book or article): Langholz, B. and Kronmal, R. A. (1991), Tests of distributional hypotheses with nuisance parameters using Fourier series,
//						   Journal of the American Statistical Association, Vol. 86, pp. 1077-1084.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat57(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$K_1$";
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
      double plaplace(double y);
      double *Y;
      Y = new double[n];
      double statK1, muhat, bhat, tmp = 0.0, tmp2 = 0.0, sumC = 0.0, sumS = 0.0, twopi;
      
      
      // calculate mu^ and b^ by using the method of moments 
      // mu^ = mean(X)
      // b^ = sqrt(sum((X-Xbar)^2)/(2n))
      
      // calculate mu^
      for (i = 0; i < n; i++) {
	tmp = tmp + x[i];
      }
      muhat = tmp / (double)n;
	
      // calculate b^
      for (i = 0; i < n; i++) {
	tmp2 = tmp2 + R_pow(x[i] - muhat, 2.0);
      }
      bhat = sqrt(tmp2 / (double)(2 * n));
      
      // generate vector Y
      twopi = 2.0 * M_PI;
      for (i = 0; i < n; i++) {
	Y[i] = twopi * plaplace((x[i] - muhat) / bhat);
      }
	
      // calculate statK1
      for (i = 0; i < n; i++) { 
	sumC = sumC + cos(Y[i]);
	sumS = sumS + sin(Y[i]);
      }
      sumC = sumC / (double)n;
      sumS = sumS / (double)n;
      
      statK1 = 2.26 * (double)n * (R_pow(sumC, 2.0) + R_pow(sumS, 2.0));
	
      statistic[0] = statK1; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue57.cpp"
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
      
    }
    
    // We return
    return;
      
  }
  
}
