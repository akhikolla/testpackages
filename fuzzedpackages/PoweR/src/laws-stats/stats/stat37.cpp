// Title: The Z_{EPD} test for normality
// Ref. (book or article): A Powerful and Interpretable Alternative to the Jarque-Bera Test of Normality Based on 2nd-Power Skewness and Kurtosis, using the Rao’s score test on the APD family, Journal of Applied Statistics

#include <R.h>
#include "Rmath.h"
#include <assert.h>

extern "C" {

  void stat37(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$Z_{EPD}$";
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
    //	double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
	double digamma(double x);	// digamma function in R extensions
    double *Y;
    Y = new double[n];   	
	double stat37, pi, euler, alphaopt, esplim, espopt, nvarlim, varopt, meanX=0.0, varX=0.0, stdX, deltahat=0.0;
	// double skewness=0.0, kurtosis=0.0;
	
	pi = 4.0*atan(1.0); 		// or use pi = M_PI, where M_PI is defined in math.h
	euler = -digamma(1.0);
	alphaopt = -0.06+2.1/R_pow((double)n,0.67);
	esplim = (R_pow(2.0-log(2.0)-euler,-0.06) - 1.0)/(-0.06);
    nvarlim = 2.0*(0.75*R_pow(pi,2.0)-7.0)/R_pow(2.0-log(2.0)-euler,2.0*(1.0-(-0.06)));
    espopt = esplim-1.32/R_pow((double)n,0.95);
	varopt = nvarlim-3.78/R_pow((double)n,0.733);
	
	// calculate meanX and varX
	for (i=0; i<n; i++) {
	  meanX = meanX + x[i];
	}
	meanX = meanX/(double)n;
	
	for (i=0; i<n; i++) {
	  varX = varX + R_pow(x[i]-meanX,2.0);
	}
	// Here we take varX = varXn for simplification
	// varX = varX/((double)(n-1)); 
	// varXn = varX*((double)(n-1))/(double)n;
	varX = varX/(double)n;
    stdX = sqrt(varX);
	
	// calculate skewness and kurtosis
	// for (i=0; i<n; i++) {
	  
	  // Y[i] = (x[i]-meanX)/stdX;
	  // skewness = skewness + R_pow(Y[i],3.0);
	  // kurtosis = kurtosis + R_pow(Y[i],4.0);
	
	// }
	// skewness = skewness/(double)n;
	// kurtosis = kurtosis/(double)n;
	
	// calculate delta^
	for (i=0; i<n; i++) {
	  Y[i] = (x[i]-meanX)/stdX;
	  if (Y[i] != 0.0) {
	    Y[i] = R_pow(Y[i],2.0);
	  } else Y[i] = 1.0;
	  deltahat = deltahat + Y[i]*log(Y[i]);
	}
	deltahat = deltahat/(double)n;
	
	// calculate stat37

	stat37 = ((R_pow(deltahat,alphaopt)-1.0)/alphaopt - espopt)/R_pow(varopt/(double)n,0.5);
	
    statistic[0] = stat37; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue37.cpp"
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
