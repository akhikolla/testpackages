// Title: Statistique de test de Chen-Shapiro
// Ref. (book or article): Chen, L. and Shapiro, S.S (1995), An alternative test for normality based on normalized spacings, 
//						   Journal of Statistical Computation and Simulation, Vol. 53,pp. 269-288.

// A VOIR !!! ON N'A PAS LA P-VALEUR!!

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat26(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
 // Here, INDICATE the name of your statistic
      const char *nom = "$CS$";
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
    double qnorm(double p, double mean, double sd, int lower_tail, int log_p);
    //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double *xs, *M;
    M = new double[n];
    xs = new double[n];
    double varX=0.0, meanX=0.0, statCS=0.0;
 

    for (i=0;i<=(n-1);i++) xs[i] = x[i];
    R_rsort (xs,n); // We sort the data
    for (i=0;i<=(n-1);i++) meanX = meanX + xs[i];
    meanX =meanX/(double)n;
    for (i=0;i<=(n-1);i++) varX = varX + R_pow(xs[i],2.0);
    varX = (varX - ((double)n)*R_pow(meanX,2.0))/(double)(n-1); 


    for (i=1;i<=n;i++) M[i-1] = qnorm(((double)i-0.375)/((double)n+0.25),0.0,1.0,1,0);

    for (i=0;i<=(n-2);i++) {
      statCS = statCS + (xs[i+1]-xs[i])/(M[i+1]-M[i]);
    }

    statCS = statCS/((double)(n-1)*sqrt(varX));

    statCS = sqrt((double)n) * (1.0 - statCS);

    statistic[0] = statCS; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue26.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
      for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values of statistic!
      } else {
	//	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }

// If applicable, we free the unused array of pointers
  delete[] M;
  delete[] xs;

}

// We return
  return;
   
        
  }
  
}
