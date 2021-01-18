// Title: Test combiné de Zhang Q-Q*
// Ref. (book or article): 
// Attention! Il s'agit d'une combinaison de 2 tests. Par conséquent, les valeurs critiques sont ici inutiles.

// A VOIR !!! ON N'A PAS LA P-VALEUR!!

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat28(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
 // Here, INDICATE the name of your statistic
      const char *nom = "Q-Q*";
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
      void stat27(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
      void stat34(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
      double pchisq(double x, double df, int lower_tail, int log_p);
      int *pvalcomp27, *pvalcomp34;
      double *pvalue27, *pvalue34, *statistic27, *statistic34;
      pvalue27 = new double[1];
      pvalue34 = new double[1];
      pvalcomp27 = new int[1];
      pvalcomp34 = new int[1];
      pvalue27[0] = 1.0;
      pvalue34[0] = 1.0;
      pvalcomp27[0] = 1;
      pvalcomp34[0] = 1;

      statistic27 = new double[1];
      statistic34 = new double[1];
      int *decision27, *decision34;
      decision27 = new int[nblevel[0]];
      decision34 = new int[nblevel[0]];
      char **name1;
      name1 = new char*[50];
      for (i=0;i<=49;i++) name1[i] =  new char[1];
      double pval1, pval2, stat;
      int *getname27, *getname34, *nbparamstat27, *nbparamstat34;
      getname27 = new int[1];
      getname27[0] = 0;
      getname34 = new int[1];
      getname34[0] = 0;
      nbparamstat27 = new int[1];
      nbparamstat27[0] = 0;
      nbparamstat34 = new int[1];
      nbparamstat34[0] = 0;
      double *paramstat27, *paramstat34;
      paramstat27 = new double[0];
      paramstat34 = new double[0];


      stat27(x,xlen,level,nblevel,name1,getname27,statistic27,pvalcomp27,pvalue27,critvalL,critvalR,usecrit,alter,decision27,paramstat27,nbparamstat27); // stat Q de Zhang

      if (pvalue27[0] > 0.5) pval1 = 1.0 - pvalue27[0]; else pval1 = pvalue27[0];

      stat34(x,xlen,level,nblevel,name1,getname34,statistic34,pvalcomp34,pvalue34,critvalL,critvalR,usecrit,alter,decision34,paramstat34,nbparamstat34); // stat Q* de Zhang

      if (pvalue34[0] > 0.5) pval2 = 1.0 - pvalue34[0]; else pval2 = pvalue34[0];

      stat = -2.0 * (log(pval1) + log(pval2)); // Combinaison des valeurs-p (Fisher, 1932)

      statistic[0] = stat; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue28.cpp"
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
      for (i=1;i<=50;i++) {
	delete[] *(name1+i-1);
      }
      delete[] name1;
      delete[] pvalue27;
      delete[] pvalue34;
      delete[] pvalcomp27;
      delete[] pvalcomp34;
      delete[] decision27;
      delete[] decision34;
      delete[] statistic27;
      delete[] statistic34;
      delete[] getname27;
      delete[] getname34;
      delete[] nbparamstat27;
      delete[] nbparamstat34;
      delete[] paramstat27;
      delete[] paramstat34;

}

// We return
    return;
   
        
  }
  
}
