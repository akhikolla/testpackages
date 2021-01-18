// Title: The Desgagné-Micheaux-Leblanc statistic for the Laplace distribution
// Ref. (book or article): Desgagné, A., Lafaye de Micheaux, P. and Leblanc, A., - unpublished document. 

#include <R.h>
#include "Rmath.h"
#include <assert.h>

extern "C" {

  void stat61(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$DLLap1$";
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
    void R_rsort(double* x, int n);
    //    double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
    double digamma(double x);	// digamma function in R extensions
    double *Y;
    Y = new double[n];   	
    double statLM, pi, euler, alphaopt, esplim, espopt, nvarlim, varopt, muhat=0.0, bhat=0.0, deltahat=0.0;
    
    pi = 4.0*atan(1.0); 		// or use pi = M_PI, where M_PI is defined in math.h
    euler = -digamma(1.0);
    alphaopt = 0.142+3.3/(double)n;
    esplim = (R_pow(1.0-euler,0.142)-1.0)/(0.142);
    nvarlim = (R_pow(pi,2.0)/3.0 - 3.0)/R_pow(1.0-euler,2.0*(1.0-(0.142)));
    espopt = esplim/(1.0+0.87/(double)n);
    varopt = nvarlim/(1.0+8.84/(double)n);
    
	// calculate mu^ and b^ by using the maximum likelihood estimators 
	// mu^ = the sample median
	// b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
	
	// calculate mu^
    R_rsort(x,n); 			// we sort the data from gensample
    if(n % 2 == 0) {		// check if n is divisible by 2
      muhat = (x[n/2-1] + x[n/2])/2.0;
    } else {
      muhat = x[n/2];
    }
    
	// calculate b^
    for (i=0;i<n;i++) {
      bhat = bhat + fabs(x[i] - muhat);
    }
    bhat = bhat/(double)n;

	
	// calculate delta^
    for (i=0; i<n; i++) {
      Y[i] = (x[i]-muhat)/bhat;
      if (Y[i] != 0.0) {
	Y[i] = fabs(Y[i]);
      } else Y[i] = 1.0;
      deltahat = deltahat + Y[i]*log(Y[i]);
    }
    deltahat = deltahat/(double)n;
	
	// calculate statLM

    statLM = ((R_pow(deltahat,alphaopt)-1.0)/alphaopt - espopt)/R_pow(varopt/(double)n,0.5);
    
    statistic[0] = statLM; // Here is the test statistic value
    
    if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue61.cpp"
    }

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
	
      if (usecrit[0] == 1) { // We use the provided critical values
	if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	} else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	} else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
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
