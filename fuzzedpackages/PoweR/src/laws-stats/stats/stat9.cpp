// Title: Statistique de test de Gel-Gastwirth robust Jarque-Bera test
// Ref. (book or article): Gel, Yulia R. and Gastwirth, Joseph L. (2008), A robust modification of the Jarque-Bera test of normality. 
//						   Economics Letters, Vol. 99, N0. 1, pp. 30-32.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat9(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$RJB$";
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
      double *xtmp;
      xtmp = new double[n];
      double m3=0.0, m4=0.0, meanX=0.0, statRJB, Jn=0.0, M;

      for (i=0;i<=(n-1);i++) meanX = meanX + x[i];
      meanX = meanX/(double)n;

      for (i=0;i<=(n-1);i++) {
	xtmp[i] = x[i];
	m3 = m3 + R_pow(x[i]-meanX,3.0);
	m4 = m4 + R_pow(x[i]-meanX,4.0);
      }
      m3 = m3/(double)n;
      m4 = m4/(double)n;
      
      R_rsort(xtmp,n);
      if ((n%2) == 0) M = (xtmp[n/2]+xtmp[n/2-1])/2.0; else M = xtmp[n/2]; // sample median
      for (i=0;i<=(n-1);i++) Jn = Jn + fabs(x[i]-M);
      Jn = sqrt(M_PI/2.0)*Jn/((double)n);
      statRJB = ((double)n)*R_pow(m3/(R_pow(Jn,3.0)),2.0)/6.0+((double)n)*R_pow(m4/R_pow(Jn,4.0)-3.0,2.0)/64.0;
      
      statistic[0] = statRJB; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue9.cpp"
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
  
}
