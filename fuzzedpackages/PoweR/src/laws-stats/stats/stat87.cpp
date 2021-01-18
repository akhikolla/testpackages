// Title: Statistique de test de Epps-Pulley with \alpha integrated out
// Ref. (book or article): Epps, T.W. and Pulley, L.B. (1983), A test of normality based on empirical characteristic function, 
//						   Biometrika, Vol. 70, No. 3, pp. 723-726.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat87(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T_n$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 0;
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
    double K = 1.0, XimXj, XimXbar;

	
    if (n>3) {
// Computation of the value of the test statistic
      double statTEP, term1=0.0, term2=0.0, xbar=0.0, S2=0.0, S, sMPI2, s1pK2;
      //double pnorm(double q, double mean, double sd, int lower_tail, int log_p);

	  // calculate xbar
	  for (i=0;i<n;i++) {
		xbar = xbar + x[i];
	  }
	  xbar = xbar/(double)n;
	  
	  // calculate S^2
	  for (i=0;i<n;i++) {
		S2 = S2 + R_pow(x[i]-xbar,2.0);
	  }
	  S2 = S2/(double)n;
	  S = sqrt(S2);
	  s1pK2 = sqrt(1 + K * K);

	  sMPI2 = sqrt(2.0 * M_PI);

	  // calculate statTEP
	  for (i=0;i<n;i++) {
	    for (j=0;j<n;j++) {
	      XimXj = fabs(x[i]-x[j]);
	      term1 = term1 + K * exp(-XimXj * XimXj / (2.0 * K * K * S2)) + sMPI2 * XimXj * (Rf_pnorm5(XimXj / (K * S), 0.0, 1.0, 1, 0) - 1.0) / S;
	    }
	  }
	  term1 = term1 / R_pow((double)n, 2.0);
	  
	  for (i=0;i<n;i++) {
	    XimXbar = x[i] - xbar;
	    term2 = term2 - exp(-XimXbar * XimXbar / (2.0 * S2)) + exp(-XimXbar * XimXbar / (2.0 * (1 + K * K) * S2)) * s1pK2 + sMPI2 * XimXbar * (Rf_pnorm5(XimXbar / (S * s1pK2), 0.0, 1.0, 1, 0) - Rf_pnorm5(XimXbar / S, 0.0, 1.0, 1, 0) ) / S;
	  }
	  term2 = 2.0 * term2 / ((double)n);
	  
	  statTEP = term1 - term2 - sqrt(2.0) + sqrt(2.0 + K * K);
	  	  
	  // end of debug
	  
      statistic[0] = ((double)n) * statTEP; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue87.cpp"
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

}

// We return
   return;
   
        
  }
  
}
