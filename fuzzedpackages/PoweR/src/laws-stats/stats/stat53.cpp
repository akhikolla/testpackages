// Title: The Choi-Kim statistic for the Laplace distribution
// Ref. (book or article): Choi, B. and Kim, K. (2006), Testing goodness-of-fit for Laplace distribution based on maximum entropy,
//						   Statistics, Vol. 40, No. 6, pp. 517-531.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat53(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 4;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{m,n}^{C}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	double mTC[] = {1.0,2.0, 2.0,3.0, 3.0,4.0, 4.0,3.0,2.0,3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
			4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
			5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	if (n > 50) {
	  paramstat[0] = 6.0;	
	} else {
	  paramstat[0] = mTC[n-4];
	}
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
	// if parstat=NULL, par = 1
	// table 4 from choi2006 for $T_{m,n}^{C}$
	double mTC[] = {1.0, 
	                2.0, 2.0, 
	                3.0, 3.0, 
				    4.0, 4.0,
					3.0,
					2.0,
					3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
					4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
				    5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	
	double par;
	if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      if (n > 50) {
	    par = 6.0;
		paramstat[0] = 6.0;	
	  } else {
	    par = mTC[n-4];
        paramstat[0] = mTC[n-4];
	  }
    } else if (nbparamstat[0] == 1) {
      par = paramstat[0];
    } else {
      return;
    }		
	

    // Here m < (n/2)
    if (n>3 && par<((double)n/2.0)) {
// Computation of the value of the test statistic
    void R_rsort(double* x, int n);
    double *Y;
    Y = new double[n];   	
    double statTCmn, tmp=0.0, bhat, muhat=0.0, sum1=0.0, sum2=0.0, prod=1.0, Zbar=0.0;
	double *Z;
	Z = new double[2*(int)par+1];
	
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
	
	// generate vector Y where the transformed data Yj = Xj - mu^, j=1,2,...,n
   	for (i=0;i<n;i++) {
	  Y[i] = x[i] - muhat;
	}
	R_rsort(Y,n);			// we sort the data, NEEDED

	// calculate statTCmn
	// remarks: Y[i-m] = Y[1] for i <= m
	//          Y[i+m] = Y[n] for i >= n-m
	// so in C, it looks like below
	
	// for (i=0; i<(n-par); i++) {
	  // if (i < par) { 
	    // Z2[i] = Y[0];
	  // } else Z2[i] = Y[i-par];
	  // if (i > (n-par-1)) {
	    // Z1[i] = Y[n-1];
	  // } else Z1[i] = Y[i+par];
	  // prod *= (Z1[i] - Z2[i]);
    // }

	// for (i=0; i<=(n-(int)par-1); i++) {
	  // Z1[i] = Y[i+(int)par];
	  // if (i <= ((int)par-1)) {
		// Z2[i] = Y[0];
	  // } else {
		// Z2[i] = Y[i-(int)par];
		// }
	  // prod *= (Z1[i] - Z2[i]);
	// }
	
	for (i=0; i<n; i++) {
	  for (j=i-(int)par; j<=(i+(int)par); j++) {
		if (j<=0) Z[j+(int)par-i] = Y[0];
	    if (j>0 && j<=(n-1)) Z[j+(int)par-i] = Y[j];
	    if (j>(n-1)) Z[j+(int)par-i] = Y[n-1];
		Zbar += Z[j+(int)par-i];	
	  }
	  Zbar = Zbar/(2.0*par+1.0);
	  for (j=i-(int)par; j<=(i+(int)par); j++) {
	    sum1 += ((((double)j+1.0)/(double)n)-(((double)i+1.0)/(double)n))*(Z[j+(int)par-i]-Zbar);
		sum2 += R_pow(Z[j+(int)par-i]-Zbar,2.0);
	  }
	  prod *= (sum1/sum2);
	  // reinitialize Zbar, sum1 and sum2
	  Zbar = 0.0;
	  sum1 = 0.0;
	  sum2 = 0.0;
	}
	
	statTCmn = R_pow(prod,-1.0/(double)n)/bhat;
	
    statistic[0] = statTCmn; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue53.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values of statistic!
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
        }
    }
    
// If applicable, we free the unused array of pointers
    delete[] Y;
	delete[] Z;

}

// We return
    return;
   
        
  }

  
}
