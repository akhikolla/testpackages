// Title: Statistique de test de Epps-Pulley
// Ref. (book or article): Epps, T.W. and Pulley, L.B. (1983), A test of normality based on empirical characteristic function, 
//						   Biometrika, Vol. 70, No. 3, pp. 723-726.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat31(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 1;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T^*(\\alpha)$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 1.0;
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
	
// Initialization of the parameters
    double alpha;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      alpha = 1.0;
      paramstat[0] = 1.0;
    } else if (nbparamstat[0] == 1) {
      alpha = paramstat[0];
    } else {
      return;
    }

	
    if (n>3) {
// Computation of the value of the test statistic
      double statTEP, term1 = 0.0, term2 = 0.0, xbar = 0.0, S2 = 0.0;

	  // calculate xbar
      for (i = 0; i < n; i++) {
	xbar = xbar + x[i];
      }
      xbar = xbar / (double)n;
      
	  // calculate S^2
      for (i = 0; i < n; i++) {
	S2 = S2 + R_pow(x[i] - xbar, 2.0);
      }
      S2 = S2 / (double)n;
	  
	  // calculate statTEP
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	  term1 = term1 + exp(-0.5 * R_pow(x[i] - x[j], 2.0) / (R_pow(alpha, 2.0) * S2));
	}
      }
      term1 = term1 / R_pow((double)n, 2.0);
	  
      for (i = 0; i < n; i++) {
	term2 = term2 + exp(-0.5 * R_pow(x[i] - xbar, 2.0) / (S2 * (1.0 + R_pow(alpha, 2.0))));
      }
      term2 = 2.0 * R_pow(1.0 + R_pow(alpha, -2.0), -0.5) * term2 / (double)n;
	  
      statTEP = term1 - term2 + R_pow(1.0 + 2.0 * R_pow(alpha, -2.0), -0.5);
	  
      statTEP = -log(((double)n) * fabs(statTEP)); 
      
      statistic[0] = statTEP; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
        #include "pvalues/pvalue31.cpp"
      }

// We take the decision to reject or not to reject the null hypothesis H0
      for (i=0;i<=(nblevel[0]-1);i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;   // less
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
