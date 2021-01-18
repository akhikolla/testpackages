// Title: The Morales statistic for uniformity
// Ref. (book or article): Morales, D., Pardo, L., Pardo, M. C. and Vajda, I. (2003), Limit laws for disparities of spacings,
//						   Journal of Nonparametric Statistics, 15(3), 325-342.
// M. A. Marhuenda, Y. Marhuenda, D. Morales, (2005), Uniformity tests under quantile categorization, 
// \emph{Kybernetes}, \bold{34}(6), 888--901.


#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat78(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$D_{n,m}(\\phi_\\lambda)$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
      paramstat[0] = 0.0;
      paramstat[1] = 2.0;
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

// Initialization of the parameters
    double lambda, m;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      lambda = 0.0;
	  m = 2.0;
      paramstat[0] = 0.0;
      paramstat[1] = 2.0;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      lambda = paramstat[0];
      m = 2.0;
      paramstat[1] = 2.0;
    } else if (nbparamstat[0] == 2) {
      lambda = paramstat[0];
      m = paramstat[1];
    } else {
      return;
    }
	

    if (n>3) {
// Computation of the value of the test statistic
    void R_rsort (double* x, int n);
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *U;
	double *G;
    U = new double[n];
	G = new double[n];
	double statDn, sumDn=0.0, nGm;
	

	// generate vector U
	for (i=0;i<n;i++) {
	  U[i] = punif(x[i],0.0,1.0,1,0);
	}
    R_rsort(U,n); // We sort the data
	
	// generate vector G
    for (i=0; i<n; i++) {
      if (i < (n-(int)m)) {
	    G[i] = U[i+(int)m] - U[i];
	  } else {
	    G[i] = 1.0 + U[i+(int)m-n] - U[i];
	  }
    }
	
	    
	// calculate statDn
      
	if (lambda == 0) {
	  for (i=0;i<n;i++) {
	    nGm = (double)n*G[i]/m;
	    if (G[i] > 0.0000000000000001) sumDn = sumDn + nGm*log(nGm);
	  }
	}
	  
	if (lambda == -1) {
	  for (i=0;i<n;i++) {
	    sumDn = sumDn - log((double)n*G[i]/m);
	  }
	}
	  
	if ((lambda != 0) && (lambda != -1)) {
	  for (i=0;i<n;i++) {
	    sumDn = sumDn + (R_pow((double)n*G[i]/m,lambda+1.0) - 1.0)/(lambda*(lambda+1.0));
	  }
	}

	statDn = sumDn;
	
	
    statistic[0] = statDn; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue78.cpp"
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
    delete[] U;
	delete[] G;

}

// We return
    return;
   
        
  }
  
  
  
}
