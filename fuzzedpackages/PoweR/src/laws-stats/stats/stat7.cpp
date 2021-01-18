// Title: The Jarque-Bera test
// Ref. (book or article): C.M. Jarque, A.K. Bera, A Test for Normality of Observations and Regression Residuals, International Statistical Review, Vol. 50, Number 2, pp. 163-172, 1987.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat7(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$JB$";
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
    double Rf_pchisq(double q, double df, int lower_tail, int log_p);
    double statJB, m2 = 0.0, m3 = 0.0, m4 = 0.0, meanX = 0.0;

    for (i = 0; i <= (n - 1); i++) {
	  meanX = meanX + x[i];
	}
    meanX = meanX / (double)n;
    
    for (i = 0; i <= (n - 1); i++) {
      m2 = m2 + R_pow(x[i] - meanX, 2.0);
      m3 = m3 + R_pow(x[i] - meanX, 3.0);
      m4 = m4 + R_pow(x[i] - meanX, 4.0);
    }
    m2 = m2 / (double)n;
    m3 = m3 / (double)n;
    m4 = m4 / (double)n;
	
    statJB = (double)n * (R_pow(m3, 2.0) / R_pow(m2, 3.0)) / 6.0 + (double)n * R_pow(m4 / R_pow(m2, 2.0) - 3.0, 2.0) / 24.0;

    statistic[0] = statJB; // Here is the test statistic value

    if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue7.cpp"
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
