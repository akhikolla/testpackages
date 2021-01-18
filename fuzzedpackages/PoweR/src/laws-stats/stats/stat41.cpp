// Title: Statistique de test de Spiegelhalter 
// Ref. (book or article): Spiegelhalter, D.J. (1977), A test for normality against symmetric alternatives, 
//						   Biometrika, Vol. 64, No. 2, pp. 415-418.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat41(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$S$";
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
    double statSp, varX=0.0, mean=0.0, max, min, u, g, sd, cn;
    double gammafn (double x);

    max = x[0];
    min = x[0];
    for (i=1;i<=(n-1);i++) {
      if (x[i] > max) max = x[i];
      if (x[i] < min) min = x[i];
    }
    for (i=0;i<=(n-1);i++) mean = mean + x[i];
    mean =mean/(double)n;
    for (i=0;i<=(n-1);i++) varX = varX + R_pow((x[i]-mean),2.0);
    varX = varX/(double)(n-1);
    sd = sqrt(varX);
    u = (max-min)/sd;
    g = 0.0;
    for (i=0;i<=(n-1);i++) g = g + fabs(x[i] - mean);
    g = g/(sd*sqrt((double)n)*sqrt((double)(n-1)));
    if (n <150) {
      cn = 0.5*R_pow(gammafn((double)(n+1)),(1.0/(double)(n-1)))/(double)n;
    } else {
      cn = R_pow(2.0*M_PI,1.0/(double)(2*(n-1)))*R_pow((n*sqrt((double)n))/M_E,1.0/(double)(n-1))/(2.0*M_E); // Stirling approximation
    }
    statSp = R_pow(R_pow(cn*u,(double)(-(n-1))) + R_pow(g,(double)(-(n-1))),(1.0/(double)(n-1)));
    
    statistic[0] = statSp; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue41.cpp"
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
