// Title: Statistique de test de Epps-Pulley modifiée avec notre poids volcano
// Ref. (book or article): Epps, T.W. and Pulley, L.B. (1983), A test of normality based on empirical characteristic function, 
//						   Biometrika, Vol. 70, No. 3, pp. 723-726.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat86(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T(\\alpha_1,\\alpha_2)$";
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
    double alpha1, alpha2;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      alpha1 = 1.0;
      paramstat[0] = 1.0;
    } else if (nbparamstat[0] == 1) {
      alpha1 = paramstat[0];
    } else {
      return;
    }

	
    if (n>3) {
// Computation of the value of the test statistic
      double statTEP=0.0, term1=0.0, term2=0.0, term3 = 0.0, term4 = 0.0, xbar=0.0, S2=0.0, alpha12, alpha22, alpha13, alpha23, XjmXk2, alpha12p1, alpha22p1, alpha12p152, alpha22p152;
      double alpha12p1S2, alpha22p1S2, tmp, XjmXbar2;
      int k, l1, l2;

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
	  
	  for (l1=24;l1<=124;l1++) {
	    alpha1 = ((double)l1) * 0.0080808;
	    
	    alpha12 = R_pow(alpha1, 2.0);
	    alpha13 = R_pow(alpha1, 3.0);
	    alpha12p1 = alpha12 + 1.0;
	    alpha12p1S2 = alpha12p1 * S2;
	    alpha12p152 = R_pow(alpha12p1, 2.5);

	    for (l2=24;l2<=124;l1++) {
	      alpha2 = ((double)l2) * 0.0080808;
		      //	    alpha2 = alpha1;

	      alpha22 = R_pow(alpha2, 2.0);
	      alpha23 = R_pow(alpha2, 3.0);
	      alpha22p1 = alpha22 + 1.0;
	      alpha22p1S2 = alpha22p1 * S2;
	      alpha22p152 = R_pow(alpha22p1, 2.5);

	  // calculate statTEP
	      XjmXbar2 = R_pow(x[0] - xbar, 2.0);
	      term3 = exp(-0.5 * XjmXbar2 / alpha12p1S2) * (alpha12p1S2 - XjmXbar2) / alpha12p152 + exp(-0.5 * XjmXbar2 / alpha22p1S2) * (alpha22p1S2 - XjmXbar2) / alpha22p152;
	      for (j=1;j<n;j++) {
		XjmXbar2 = R_pow(x[j] - xbar, 2.0);
		for (k=0;k<j;k++) {
		  XjmXk2 = R_pow(x[j] - x[k], 2.0);
		  term1 = term1 + exp(-0.5 * XjmXk2  / (alpha22 * S2)) * (alpha22 * S2 - XjmXk2);
		  term2 = term2 + exp(-0.5 * XjmXk2  / (alpha12 * S2)) * (alpha12 * S2 - XjmXk2);
		}
		term3 = term3 + exp(-0.5 * XjmXbar2 / alpha12p1S2) * (alpha12p1S2 - XjmXbar2) / alpha12p152 + exp(-0.5 * XjmXbar2 / alpha22p1S2) * (alpha22p1S2 - XjmXbar2) / alpha22p152;
	      }
	      term1 = 2.0 * term1 / R_pow((double)n, 2.0);
	      term2 = 2.0 * term2 / R_pow((double)n, 2.0);
	      term1 = R_pow(alpha1 / alpha2, 2.0) * (alpha1 / (S2 * (alpha13 + alpha23))) * term1;
	      term2 = R_pow(alpha2 / alpha1, 2.0) * (alpha2 / (S2 * (alpha13 + alpha23))) * term1;


	      term3 = 2.0 * alpha13 * alpha23 * term3 / S2;
	      term3 = term3 / (alpha13 + alpha23);
	      term3 = term3 / ((double)n);
	      
	      term4 = (alpha13 * alpha23 * (sqrt(alpha12 + 2.0) * alpha12 + 2.0 * sqrt(alpha12 + 2.0) + sqrt(alpha22 + 2.0) * alpha22 + 2.0 * sqrt(alpha22 + 2.0))) / (R_pow(alpha12 + 2.0, 1.5) * R_pow(alpha22 + 2.0, 1.5) * (alpha13 + alpha23));

	      tmp = (1.0 / ((double)n)) + term1 + term2 - term3 + term4;

	      //   if (tmp > statTEP2) statTEP2 = tmp;
	      statTEP = statTEP + tmp;
		    }
	  }
	  
	  // end of debug
	  
      statistic[0] = ((double)n) * statTEP ; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue86.cpp"
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
