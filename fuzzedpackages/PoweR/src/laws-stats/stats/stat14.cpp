// Title: Statistique de test de Bontemps-Medahhi BM_{3-4}
// Ref. (book or article): Bontemps, C. and Meddahi, N. (2005), \emph{Testing Normality: A GMM Approach}, 
//						   Journal of Econometrics, Vol. 124, pp. 149-186.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat14(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$BM_{3-4}$";
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
      void R_rsort (double* x, int n);
      double pchisq(double q, double df, int lower_tail, int log_p);
      double *z;
      z = new double[n];
      double varX=0.0, meanX=0.0, sdX, statBM34=0.0, tmp3=0.0, tmp4=0.0;
      
      for (i=0;i<=(n-1);i++) meanX = meanX + x[i];
      meanX =meanX/(double)n;
      for (i=0;i<=(n-1);i++) varX = varX + R_pow(x[i],2.0);
      varX = ((double)n)*(varX/(double)n - R_pow(meanX,2.0))/(double)(n-1); 
      sdX = sqrt(varX);

      for (i=0;i<=(n-1);i++) z[i] = (x[i]-meanX)/sdX;
      
      for (i=0;i<=(n-1);i++) {
		tmp3 = tmp3 + (R_pow(z[i],3.0) - 3.0*z[i])/sqrt(6.0);
		tmp4 = tmp4 + (R_pow(z[i],4.0) - 6.0*R_pow(z[i],2.0) + 3.0)/(2.0*sqrt(6.0));
      }

      statBM34 = (R_pow(tmp3,2.0) + R_pow(tmp4,2.0))/(double)n;  
      
      statistic[0] = statBM34; // Here is the test statistic value
      
if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue14.cpp"
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
      delete[] z;
    }

// We return
    return;
      
  }

   
}
