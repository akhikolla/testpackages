// Title: Statistique de test de Bonett-Seier
// Ref. (book or article): Bonett, D.G. and Seier, E. (2002), A test of normality with high uniform power, 
//						   Computational Statistics & Data Analysis, Vol. 40, pp. 435-445.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat17(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_w$";
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
      double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
      double statTw, m2=0.0, meanX=0.0, omega, term=0.0;

      for (i=0;i<=(n-1);i++) {
	    meanX = meanX + x[i];
	  }
      meanX = meanX/(double)n;
      
	  for (i=0;i<=(n-1);i++) {
		m2 = m2 + R_pow(x[i]-meanX,2.0);
		term = term + fabs(x[i]-meanX);
      }
      m2 = m2/(double)n;
      term = term/(double)n;
      
      omega = 13.29*(log(sqrt(m2))-log(term));
      statTw = sqrt((double)(n+2))*(omega-3.0)/3.54;
      
      statistic[0] = statTw; // Here is the test statistic value

	  
if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue17.cpp"
}
      
// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
	
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	  } else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	  } else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
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
