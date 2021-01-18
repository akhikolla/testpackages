// Title: The K2 test D'Agostino-Pearson
// Ref. (book or article): Warning: different from the agostino.test() function in the moments R package

// A VOIR !!! ON N'A PAS LA P-VALEUR!!

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat6(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
 // Here, INDICATE the name of your statistic
      const char *nom = "$K^2$";
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
      for (i = j; i < 50;i++) name[i][0] = space[0];
      return;
    }
  
    
    if (n > 3) {
// Computation of the value of the test statistic
      double pchisq(double q, double df, int lower_tail, int log_p);
      double varX = 0.0, meanX = 0.0, statK2 = 0.0, skew = 0.0, kurt = 0.0, tmp, yy, y, beta2, w2, c, term1, eb2, vb2, srbeta1, A, term2;
 
      for (i = 0; i < n; i++) meanX = meanX + x[i];
      meanX = meanX / (double)n;
      for (i = 0; i < n; i++) varX = varX + R_pow(x[i], 2.0);
      varX = varX / (double)n -R_pow(meanX, 2.0); // !! var. pop. ici
      for (i = 0; i < n; i++) {
	  tmp = x[i] - meanX;
	  skew = skew + R_pow(tmp, 3.0);   
	  kurt = kurt + R_pow(tmp, 4.0);
      }
      skew = skew / (R_pow(varX, 3.0 / 2.0) * (double)n);
      kurt = kurt / (R_pow(varX, 2.0) * (double)n);

      yy = skew * sqrt(((double)(n + 1) * (n + 3)) / ((double)(6 * (n - 2))));
      //      beta2 = ((double)(3*(n*n+27*n-70)*(n+1)*(n+3)))/((double)((n-2)*(n+5)*(n+7)*(n+9))); // !! On débordait la capacité des entiers 2^32 ici !!
      beta2 = (3.0 * (double)(n * n + 27 * n - 70) * (double)(n + 1) * (double)(n + 3)) / ((double)(n - 2) * (double)(n + 5) * (double)(n + 7) * (double)(n + 9));
      w2 = -1.0 + sqrt(2.0 * (beta2 - 1.0));
      c = sqrt(2.0 / (w2 - 1.0));
      term1 = (log(yy / c + sqrt(R_pow((yy / c), 2.0) + 1.0))) / sqrt(log(sqrt(w2)));

      eb2 = ((double)(3 * (n - 1))) / ((double)(n + 1));
      vb2 = ((24 * (double)n * (double)(n - 2) * (double)(n - 3))) / ((double)(n + 1) * (double)(n + 1) * (double)(n + 3) * (double)(n + 5));
      y = (kurt - eb2) / sqrt(vb2);
      srbeta1 = ((6.0 * (double)(n * n - 5 * n + 2)) / ((double)(n + 7) * (double)(n + 9))) * sqrt((6.0 * (double)(n + 3) * (double)(n + 5)) / ((double)n * (double)(n - 2) * (double)(n - 3)));
      A = 6.0 + (8.0 / srbeta1) * (2.0 / srbeta1 + sqrt(1.0 + 4.0 / (srbeta1 * srbeta1)));
      term2 = ((1.0 - 2.0 / (9.0 * A)) - R_pow((1.0 - 2.0 / A) / (1.0 + y * sqrt(2.0 / (A - 4.0))), (1.0 / 3.0))) / sqrt(2.0 / (9.0 * A));
      statK2 = term1 * term1 + term2 * term2;
      
      
      statistic[0] = statK2; // Here is the test statistic value
      

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue6.cpp"
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
