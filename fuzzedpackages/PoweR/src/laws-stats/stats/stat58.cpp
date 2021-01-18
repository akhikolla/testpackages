// Title: The Kundu statistic for the Laplace distribution
// Ref. (book or article): Kundu, Debasis (2005), Discriminating between Normal and Laplace distributions,
//						   Advances in ranking and selection, multiple comparisons, and reliability, 65--79, 
//						   Stat. Ind. Technol., Birkhäuser Boston, Boston, MA. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat58(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T$";
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
    // double qlaplace(double p, double location, double scale);
	double qnorm(double p, double mu, double sigma, int lower_tail, int log_p);
    double statT, etahat=0.0, thetahat, muhat=0.0, sigmahat, tmp=0.0, tmp2=0.0, tmp3=0.0, pi;
	
	pi = 4.0*atan(1.0); 		// or use pi = M_PI, where M_PI is defined in math.h

	// calculate eta^ and theta^ by using the maximum likelihood estimators 
	// eta^ = the sample median
	// theta^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
	
	// calculate eta^
	R_rsort(x,n); 			// we sort the data from gensample
	if(n % 2 == 0) {		// check if n is divisible by 2
	  etahat = (x[n/2-1] + x[n/2])/2.0;
    } else {
      etahat = x[n/2];
	}
	
	// calculate theta^
	for (i=0;i<n;i++) {
	  tmp = tmp + fabs(x[i] - etahat);
	}
	thetahat = tmp/(double)n;
	
	// calculate mu^
	for (i=0;i<=(n-1);i++) {
	  tmp2 = tmp2 + x[i];
    }
	muhat = tmp2/(double)n;
	
	// calculate sigma^	
	for (i=0; i<n; i++) {
	  tmp3 = tmp3 + R_pow(x[i]-muhat,2.0);
	}
    sigmahat = sqrt(tmp3/(double)n);	
	
	
	// calculate statT
		
	// statT = (double)n*log(2.0)/2.0 - (double)n*log(pi)/2.0 + (double)n*log(thetahat) - (double)n*log(sigmahat) + (double)n/2.0;
	statT = (double)n*(log(2.0)/2.0 - log(pi)/2.0 + log(thetahat) - log(sigmahat) + 1.0/2.0);
	
    statistic[0] = statT; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue58.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    // for (i=0;i<=(nblevel[0]-1);i++) {
	  
	  // if T > 0 => Normal
	  // if T <= 0 => Laplace
	  // if (statistic[0] > -(double)n*0.0723649+qnorm(1.0-level[i],0.0,1.0,1,0)*sqrt((double)n*0.25)) decision[i] = 1; else decision[i] = 0; 
   
    // }
	
	for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	    if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		  // If critval = NULL, we use the expression from Kundu (2005) pp. 72 - Test 2
		  if (statistic[0] > -(double)n*0.0723649+qnorm(1.0-level[i],0.0,1.0,1,0)*sqrt((double)n*0.25)) decision[i] = 1; else decision[i] = 0;  
        }
    }
    
// If applicable, we free the unused array of pointers

}

// We return
    return;
   
        
  }
  
  // The quantile function for the Laplace distribution 
  // double qlaplace(double p, double location, double scale) {
    
	// return (p < 0.5) ? location + scale*log(2.0*p) : location - scale*log(2.0*(1.0-p));
  
  // }
  
  
  
}
