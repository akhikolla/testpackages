// Title: Statistique de test de Hosking based on trimmed L-moments T_{Lmom}^{(2)}
// Ref. (book or article): Hosking, J.R.M. (1990). "L-moments: analysis and estimation of distributions using linear combinations of order statistics",
//						   Journal of the Royal Statistical Society, Series B 52: 105-124.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat12(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{Lmom}^{(2)}$";
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
      double pstarmod2(int r, int k, int i);
      double *xtmp;
      xtmp = new double[n];
      double statTLmom2, l22=0.0, l32=0.0, l42=0.0, tau32, tau42;
      double mutau42, vtau32, vtau42;
      

      for (i=0;i<=(n-1);i++) xtmp[i] = x[i];
      
      R_rsort(xtmp,n); // We sort the data

      for (i=2;i<=(n-1);i++) {

	l22 = l22 + xtmp[i-1]*pstarmod2(2,n,i);
	l32 = l32 + xtmp[i-1]*pstarmod2(3,n,i);
	l42 = l42 + xtmp[i-1]*pstarmod2(4,n,i);

      }

      l22 = l22/(2.0*choose((double)n,6.0));
      l32 = l32/(3.0*choose((double)n,7.0));
      l42 = l42/(4.0*choose((double)n,8.0));

      tau32 = l32/l22;
      tau42 = l42/l22;

      if (1<=n && n<=25) {
	mutau42 = 0.044174;
	vtau32 = 0.0086570;
	vtau42 = 0.0042066;
      }

      if (25<n && n<=50) {
	mutau42 = 0.040389;
	vtau32 = 0.0033818;
	vtau42 = 0.0013301;
      }
 
      if (50<n) {
	mutau42 = 0.039030;
	vtau32 = 0.0015120;
	vtau42 = 0.00054207;
      }

      statTLmom2 = R_pow(tau32,2.0)/vtau32 + R_pow(tau42 - mutau42,2.0)/vtau42;
      
      statistic[0] = statTLmom2; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue12.cpp"
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


  double pstarmod2(int r, int n, int i) {

    double choose(double n, double k);
    double res=0.0;
    int k;

    for (k=0;k<=(r-1);k++) {

      res = res + R_pow(-1.0,(double)k) * choose((double)(r-1),(double)k) * choose((double)(i-1),(double)(r+2-1-k)) * choose((double)(n-i),(double)(2+k));

    }

    return(res);

  }

  
}
