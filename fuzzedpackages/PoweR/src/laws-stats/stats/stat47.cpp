// Title: The Meintanis statistic for the Laplace distribution - MO estimators
// Ref. (book or article): Meintanis, S.G. (2004), A Class of Omnibus Tests for the Laplace Distribution Based on the Empirical Characteristic Function,
//						   Communications in Statistics - Theory and Methods, Vol. 33, No. 4, pp. 925-948.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat47(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{n,a}^{(1)}-MO$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 2.0;
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
    double a;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      a = 2.0;
      paramstat[0] = 2.0;
    } else if (nbparamstat[0] == 1) {
      a = paramstat[0];
    } else {
      return;
    }
	
	
    if (n > 3) {
// Computation of the value of the test statistic
    double *Y;
    Y = new double[n];
    double statT1na, tmp = 0.0, tmp2 = 0.0, bhat, muhat, sumM1 = 0.0, sumM2 = 0.0;
    double a2 = a * a, a4 = a2 * a2, Yi2, YimYj2;
    
    // calculate mu^ and b^ by using the moments estimators 
    // mu^ = Xbar (sample mean)
    // b^ = S / sqrt(2) where S^2 = sample variance = n^{-1} \sum_{j=1}^n(X_j-\bar{X})^2
	
    // calculate mu^
    for (i = 0; i < n; i++) {
      tmp = tmp + x[i];
    }
    muhat = tmp / (double)n;
	
    // calculate b^
    for (i = 0; i < n; i++) {
      tmp2 = tmp2 + R_pow(x[i] - muhat, 2.0);
    }
    bhat = sqrt(tmp2 / (double)(2 *n));	
	
    // generate vector Y where the transformed data Yj = (Xj - mu^) / b^, j = 1, 2, ..., n
    for (i = 0; i < n; i++) {
      Y[i] = (x[i] - muhat) / bhat;
    }
	
    // calculate statT1na
    for (i = 0; i < n; i++) {
      Yi2 = R_pow(Y[i], 2.0);
      sumM1 = sumM1 + 1.0 / (a2 + Yi2) + 2.0 * (a2 - 3.0 * Yi2) / R_pow(a2 + Yi2, 3.0);
      for (j = 0; j < n; j++) {
	YimYj2 = R_pow(Y[i] - Y[j], 2.0);
	sumM2 = sumM2
	  + 1.0 / (a2 + YimYj2)
	  + 24.0 * (a4 + 5.0 * YimYj2 * YimYj2 - 10.0 * a2 * YimYj2) / R_pow(a2 + YimYj2, 5.0) 
	  + 4.0 * (a2 - 3.0 * YimYj2) / R_pow(a2 + YimYj2, 3.0);
      }
    }
    
    statT1na = (double)(2 * n) / a - 4.0 * a * sumM1 + 2.0 * a * sumM2 / (double)n;
	
    statistic[0] = statT1na; // Here is the test statistic value

    if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue47.cpp"
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
    delete[] Y;
    
    }
    
    // We return
    return;
    
  }
  
}
