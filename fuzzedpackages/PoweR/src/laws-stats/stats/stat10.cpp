// Title: Statistique de test de Hosking L-moments
// Ref. (book or article): Hosking, J.R.M. (1990). "L-moments: analysis and estimation of distributions using linear combinations of order statistics",
//						   Journal of the Royal Statistical Society, Series B 52: 105-124.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat10(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{Lmom}$";
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
      double pchisq(double q, double df, int lower_tail, int log_p);
      double *xtmp;
      xtmp = new double[n];
      double statTLmom, l2, l3, l4, tau3, tau4, b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
      double mutau4, vtau3, vtau4;
      int tmp1, tmp2, tmp3;

      for (i = 0; i < n; i++) xtmp[i] = x[i];
      
      R_rsort(xtmp, n); // We sort the data

      tmp1 = n * (n - 1);
      tmp2 = tmp1 * (n - 2);
      tmp3 = tmp2 * (n - 3);
      
      b0 = b0 + xtmp[0] + xtmp[1] + xtmp[2];
      b1 = b1 + 1.0 * xtmp[1] + 2.0 * xtmp[2];
      b2 = b2 + 2.0 * xtmp[2];
      
      for (i = 3; i < n; i++) {
	b0 = b0 + xtmp[i];
	b1 = b1 + ((double)i) * xtmp[i];
	b2 = b2 + (double)(i * (i - 1)) * xtmp[i];
	b3 = b3 + (double)(i * (i - 1) * (i - 2)) * xtmp[i];
      }

      b0 = b0 / (double)n;
      b1 = b1 / (double)tmp1;
      b2 = b2 / (double)tmp2;
      b3 = b3 / (double)tmp3;

      l2 = 2.0 * b1 - b0;  // l2 = p_{1,0}b0 + p_{1,1}b1
      l3 = 6.0 * b2 - 6.0 * b1 + b0; // l3 = p_{2,0}b0 + p_{2,1}b1 + p_{2,2}b2
      l4 = 20.0 * b3 - 30.0 * b2 + 12.0 * b1 - b0; // l4 = p_{3,0}b0 + p_{3,1}b1 + p_{3,2}b2 + p_{3,3}b3

      tau3 = l3 / l2;
      tau4 = l4 / l2;

      // Table 1 page 563 of Romao et al.
      if ((1 <= n) && (n <= 25)) {
	mutau4 = 0.12383;
	vtau3 = 0.0088038;
	vtau4 = 0.0049295;
      }

      if ((25 < n) && (n <= 50)) {
	mutau4 = 0.12321;
	vtau3 = 0.0040493;
	vtau4 = 0.0020802;
      }
      
      if (n > 50) {
	mutau4 = 0.12291;
	vtau3 = 0.0019434;
	vtau4 = 0.00095785;
      }

      // I corrected the error in Romao et al. (missing squares)
      statTLmom = R_pow(tau3, 2.0) / vtau3 + R_pow(tau4 - mutau4, 2.0) / vtau4;
      
      statistic[0] = statTLmom; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue10.cpp"
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
      delete[] xtmp;
    }
    
// We return
    return;
       
  }
 
}
