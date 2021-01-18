// Title: More Light on the Kurtosis and Related Statistics (for the Laplace distribution)
// Ref. (book or article): Hogg, R. V. 1972. More Light on the Kurtosis and Related Statistics. Journal
//                         of the American Statistical Association, 67(338), 422-424.
#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat93(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 0;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$Ho_K$";
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
    
    // Initialization of the parameters
	
	
    if (n > 3) {
      // Computation of the value of the test statistic
      double xbar = 0.0, S2 = 0.0, statHoK, temp, num = 0.0;
      // calculate xbar
      for (i = 0; i < n; i++) {
	xbar = xbar + x[i];
      }
      xbar = xbar / (double)n;
      
      // calculate S^2
      for (i = 0; i < n; i++) {
	temp = R_pow(x[i] - xbar, 2.0);
	S2 = S2 + temp;
	num = num + R_pow(temp, 2.0);
      }

      statHoK = (double)n * num / R_pow(S2, 2.0);

      statistic[0] = statHoK; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue93.cpp"
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
