// Title: The Doornik-Hansen test
// Ref. (book or article): J.A. Doornik, H. Hansen, An Omnibus Test for Univariate and Multivariate Normality, 1994.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat8(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$DH$";

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
      double Rf_pchisq(double q, double df, int lower_tail, int log_p);
      double statDH, m2 = 0.0, m3 = 0.0, m4 = 0.0, meanX = 0.0, skew, kurt, yy, beta2, w2, c, term1, z2, xi, d, a, k;
      int n2 = n * n;

      for (i = 0; i< n; i++) meanX = meanX + x[i];
      meanX = meanX / (double)n;
      for (i = 0; i < n; i++) {
		m2 = m2 + R_pow(x[i] - meanX, 2.0);
		m3 = m3 + R_pow(x[i] - meanX, 3.0);
		m4 = m4 + R_pow(x[i] - meanX, 4.0);
      }
      m2 = m2 / (double)n;
      m3 = m3 / (double)n;
      m4 = m4 / (double)n;
      
      skew = m3 / R_pow(m2, 3.0 / 2.0);
      kurt = m4 / R_pow(m2, 2.0);
      
      yy = skew * sqrt(((double)(n + 1) * (double)(n + 3)) / ((double)(6 * (n - 2))));
      beta2 = (3.0 * ((double)n2 + 27 * (double)n - 70.0) * (double)(n + 1) * (double)(n + 3)) / ((double)(n - 2) * (double)(n + 5) * (double)(n + 7) * (double)(n + 9));
      w2 = -1.0 + sqrt(2.0 * (beta2 - 1.0));
      c = sqrt(2.0 / (w2 - 1.0));
      term1 = (log(yy / c + sqrt(R_pow((yy / c), 2.0) + 1.0))) / sqrt(log(sqrt(w2)));
      
      d = (double)(n - 3) * (double)(n + 1) * (double)(n2 + 15 * n - 4);
      a = ((double)(n - 2) * (double)(n + 5) * (double)(n + 7) * (double)(n2 + 27 * n - 70)) / (6.0 * d);
      c = ((double)(n - 7) * (double)(n + 5) * (double)(n + 7) * (double)(n2 + 2 * n - 5)) / (6.0 * d);
      k = ((double)(n + 5) * (double)(n + 7) * (double)(n * n2 + 37 * n2 + 11 * n - 313)) / (12.0 * d);  
      a = a + skew * skew * c;
      xi = (kurt - 1.0 - skew * skew) * 2.0 * k;
      xi = fabs(xi);

      z2 = (R_pow(xi / (2.0 * a), 1.0 / 3.0) - 1.0 + 1.0 / (9.0 * a)) * sqrt(9.0 * a);
      
      statDH = term1 * term1 + z2 * z2;
      
      statistic[0] = statDH; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue8.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i < nblevel[0];i++) {
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
