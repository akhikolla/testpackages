// Title: Test combiné de Bonett-Seier T_w et de Brys-Hubert-Struyf MC-LR
// Attention! Il s'agit d'une combinaison de 2 tests. Par conséquent, les valeurs critiques sont ici inutiles.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat18(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{MC-LR}-T_w$";
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
      void stat16(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
      void stat17(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
      double pchisq(double x, double df, int lower_tail, int log_p);
      double *pvalue16, *pvalue17, *statistic16, *statistic17;
      int *pvalcomp16, *pvalcomp17;
      pvalue16 = new double[1];
      pvalue16[0] = 0.0;
      pvalue17 = new double[1];
      pvalue17[0] = 0.0;
      pvalcomp16 = new int[1];
      pvalcomp17 = new int[1];
      pvalcomp16[0] = 1;
      pvalcomp17[0] = 1;
      statistic16 = new double[1];
      statistic17 = new double[1];
      int *alter17, *decision16, *decision17;
      alter17 = new int[1];
      alter17[0] = 0;
      decision16 = new int[nblevel[0]];
      decision17 = new int[nblevel[0]];
      char **name1;
      name1 = new char*[50];
      for (i=0;i<=49;i++) name1[i] =  new char[1];
      double pval1, pval2, stat;
      int *getname16, *getname17, *nbparamstat16, *nbparamstat17;
      getname16 = new int[1];
      getname16[0] = 0;
      getname17 = new int[1];
      getname17[0] = 0;
      nbparamstat16 = new int[1];
      nbparamstat16[0] = 0;
      nbparamstat17 = new int[1];
      nbparamstat17[0] = 0;
      double *paramstat16, *paramstat17;
      paramstat16 = new double[0];
      paramstat17 = new double[0];


      stat16(x,xlen,level,nblevel,name1,getname16,statistic16,pvalcomp16,pvalue16,critvalL,critvalR,usecrit,alter,decision16,paramstat16,nbparamstat16); // stat T_{MC-LR}

      if (pvalue16[0] > 0.5) pval1 = 1.0 - pvalue16[0]; else pval1 = pvalue16[0];

      stat17(x,xlen,level,nblevel,name1,getname17,statistic17,pvalcomp17,pvalue17,critvalL,critvalR,usecrit,alter17,decision17,paramstat17,nbparamstat17); // stat T_w

      if (pvalue17[0] > 0.5) pval2 = 1.0 - pvalue17[0]; else pval2 = pvalue17[0];

      stat = -2.0 * (log(pval1) + log(pval2)); // Combinaison des valeurs-p (Fisher, 1932)

      statistic[0] = stat; // Here is the test statistic value

	  
if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue18.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
      for (i=0;i<=(nblevel[0]-1);i++) {
	
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0;   // greater
	} else {
	  if (pvalue[0] < 0.5*level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
	}
      }

 // If applicable, we free the unused array of pointers
      for (i=1;i<=50;i++) {
	delete[] *(name1+i-1);
      }
      delete[] name1;
      delete[] alter17;
      delete[] pvalue16;
      delete[] pvalue17;
      delete[] pvalcomp16;
      delete[] pvalcomp17;
      delete[] decision16;
      delete[] decision17;
      delete[] statistic16;
      delete[] statistic17;
      delete[] getname16;
      delete[] getname17;
      delete[] nbparamstat16;
      delete[] nbparamstat17;
      delete[] paramstat16;
      delete[] paramstat17;


}

// We return
    return;
   
        
  }
  
}
