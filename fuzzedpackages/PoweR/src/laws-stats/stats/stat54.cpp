// Title: The Desgagné-Micheaux-Leblanc statistic for the Laplace distribution
// Ref. (book or article): Desgagné, A., Lafaye de Micheaux, P. and Leblanc, A., - unpublished document. 

#include <R.h>
#include "Rmath.h"
#include <assert.h>

extern "C" {

  void stat54(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$\\hat{G}_n$";
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
	double pchisq(double q, double df, int lower_tail, int log_p);
	double digamma(double x);	// digamma function in R extensions
    double *Y;
    Y = new double[n];   	
    double statGn, tmp=0.0, bhat, muhat=0.0, sum=0.0, pi, K, gamma, gn;
	
	
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
	  tmp = tmp + fabs(x[i] - muhat);
	}
	bhat = tmp/(double)n;
	
	// generate vector Y where the transformed data Yj = (Xj - mu^)/b^, j=1,2,...,n
   	for (i=0;i<n;i++) {
	  Y[i] = (x[i] - muhat)/bhat;
	}
	// R_rsort(Y,n);			// we sort the data, NEEDED

	// calculate statGn
	pi = 4.0*atan(1.0); 		// or use pi = M_PI, where M_PI is defined in math.h
	K = R_pow(pi,2.0)/3.0 - 3.0;
	gamma = 1.0 - digamma(2.0); // the Euler-Mascheroni constant
	
	for (i=0; i<n; i++) {
	
	  if (Y[i] == 0.0) Y[i] = 1.0; 
	  
	  sum = sum + fabs(Y[i])*log(fabs(Y[i]));
	  
	}
	gn = 1.0 - gamma - sum/(double)n;
	
	statGn = (double)n*R_pow(gn,2.0)/K;
	
    statistic[0] = statGn; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue54.cpp"
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
    delete[] Y;

}

// We return
    return;
   
        
  }

  
}
