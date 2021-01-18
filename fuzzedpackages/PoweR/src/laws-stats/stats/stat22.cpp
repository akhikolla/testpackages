// Title: Statistique de test de Shapiro-Francia 
// Ref. (book or article): (Voir le package nortest) 
//						   Shapiro, S.S. and Francia, R. (1972), An approximation analysis of variance test for normality, 
//						   Journal of the American Statistical Association, Vol. 67, pp. 215-216.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat22(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 4;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$W'$";
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
    //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double qnorm(double p, double mean, double sd, int lower_tail, int log_p);
    double *y, *xs;
    xs = new double[n];
    y = new double[n];
    double varX=0.0, varY=0.0, meanX=0.0, meanY=0.0, a, statWSF=0.0, u, v, mu, sig, z;
 
    for (i=0;i<=(n-1);i++) xs[i] = x[i];
    R_rsort (xs,n); // We sort the data
    a = 3.0/8.0;
    for (i=1;i<=n;i++) y[i-1] = qnorm(((double)i-a)/((double)n+1-2*a),0.0,1.0,1,0);
    for (i=0;i<=(n-1);i++) {
      meanX = meanX + xs[i];
      meanY = meanY + y[i];
    }
    meanX =meanX/(double)n;
    meanY =meanY/(double)n;
    for (i=0;i<=(n-1);i++) {
      varX = varX + R_pow(xs[i],2.0);
      varY = varY + R_pow(y[i],2.0);
    }
    varX = varX/(double)(n-1) -((double)n)*R_pow(meanX,2.0)/(double)(n-1);
    varY = varY/(double)(n-1) -((double)n)*R_pow(meanY,2.0)/(double)(n-1);
    for (i=0;i<=(n-1);i++) {
     statWSF = statWSF + xs[i]*y[i];
   }
   statWSF = R_pow(statWSF/(double)(n-1) - ((double)n)*meanX*meanY/(double)(n-1),2.0)/(varX*varY);
   u = log((double)n);
   v = log(u);
   mu = -1.2725 + 1.0521 * (v - u);
   sig = 1.0308 - 0.26758 * (v + 2.0/u);
   z = (log(1.0 - statWSF) - mu)/sig;

    statistic[0] = statWSF; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue22.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
   for (i=0;i<=(nblevel[0]-1);i++) {   
      if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values
      } else {
	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
   
// If applicable, we free the unused array of pointers
   delete[] xs;
   delete[] y;

}

// We return
   return;
   
        
  }
  
}
