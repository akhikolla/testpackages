// Title: The Cramér-von Mises statistic for uniformity
// Ref. (book or article): Anderson, T.W. and Darling, D.A. (1954) 'A test of goodness-of-fit', 
//                         Journal of the American Statistical Association, 49, 765-769.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat64(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$W_{n}^{2}$";
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
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *CvM;
    CvM = new double[n];
    double statCvM, sumCvM=0.0;
	
	// generate vector CvM
    for (i=0;i<n;i++) {
	  CvM[i] = punif(x[i],0.0,1.0,1,0);
	}
    R_rsort(CvM,n); // We sort the data, NO NEED SINCE plaplace(x[i]) < plaplace(x[i+1])
	
	// calculate statCvM
    for (i=1;i<=n;i++) {
	  sumCvM = sumCvM + R_pow((CvM[i-1] - (2.0*(double)i-1)/(2.0*(double)n)),2.0);
	}
    statCvM = 1.0/(12.0*(double)n) + sumCvM;

    statistic[0] = statCvM; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue64.cpp"
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
    delete[] CvM;

}

// We return
    return;
   
        
  }

  
}
