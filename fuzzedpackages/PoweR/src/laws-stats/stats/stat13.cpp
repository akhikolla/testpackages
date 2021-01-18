// Title: Statistique de test de Hosking based on trimmed L-moments T_{Lmom}^{(3)}
// Ref. (book or article): Hosking, J.R.M. (1990). "L-moments: analysis and estimation of distributions using linear combinations of order statistics",
//						   Journal of the Royal Statistical Society, Series B 52: 105-124.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat13(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{Lmom}^{(3)}$";
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
      double pstarmod3(int r, int k, int i);
      double *xtmp;
      xtmp = new double[n];
      double statTLmom3, l23=0.0, l33=0.0, l43=0.0, tau33, tau43;
      double mutau43, vtau33, vtau43;
      

      for (i=0;i<=(n-1);i++) xtmp[i] = x[i];
      
      R_rsort(xtmp,n); // We sort the data

      for (i=2;i<=(n-1);i++) {

	l23 = l23 + xtmp[i-1]*pstarmod3(2,n,i);
	l33 = l33 + xtmp[i-1]*pstarmod3(3,n,i);
	l43 = l43 + xtmp[i-1]*pstarmod3(4,n,i);

      }

      l23 = l23/(2.0*choose((double)n,8.0));
      l33 = l33/(3.0*choose((double)n,9.0));
      l43 = l43/(4.0*choose((double)n,10.0));

      tau33 = l33/l23;
      tau43 = l43/l23;

      if (1<=n && n<=25) {
	mutau43 = 0.033180;
	vtau33 = 0.0095765;
	vtau43 = 0.0044609;
      }

      if (25<n && n<=50) {
	mutau43 = 0.028224;
	vtau33 = 0.0033813;
	vtau43 = 0.0011823;
      }
 
      if (50<n) {
	mutau43 = 0.026645;
	vtau33 = 0.0014547;
	vtau43 = 0.00045107;
      }

      statTLmom3 = R_pow(tau33,2.0)/vtau33 + R_pow(tau43 - mutau43,2.0)/vtau43;
      
      statistic[0] = statTLmom3; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue13.cpp"
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
      delete[] xtmp;
    }
    
// We return
    return;
   
    
  }


  double pstarmod3(int r, int n, int i) {

    double choose(double n, double k);
    double res=0.0;
    int k;

    for (k=0;k<=(r-1);k++) {

      res = res + R_pow(-1.0,(double)k) * choose((double)(r-1),(double)k) * choose((double)(i-1),(double)(r+3-1-k)) * choose((double)(n-i),(double)(3+k));

    }

    return(res);

  }

  
}
