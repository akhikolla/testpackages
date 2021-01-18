// Title: GIVE A SHORT DESCRIPTION OF YOUR TEST.
// Ref.: GIVE A REFERENCE (BOOK OR ARTICLE FOR YOUR TEST).

#include <R.h>
#include "Rmath.h"

extern "C" {

  void statj(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}"); // MODIFY (OR REMOVE) THIS LINE TO SUIT YOUR NEEDS.

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// HERE, INDICATE THE NAME OF YOUR STATISTIC.
      const char *nom = "name_of_your_statistic";
// HERE, INDICATE THE NUMBER OF PARAMETERS OF YOUR STATISTIC (E.G., 1). WRITE 0 IF YOUR STATISTIC INVOLVES NO PARAMETER.
      nbparamstat[0] = 1;
      if (name[0][0] == '1') { // DO NOT MODIFY THIS LINE.
	// REMOVE THE LINE BELOW IF NO PARAMETER IS INVOLVED IN YOUR STATISTIC.
	// MODIFY THE LINE BELOW TO SUIT YOUR NEEDS IF ONLY ONE PARAMETER IS INVOLVED IN YOUR STATISTIC (PROVIDE THE DEFAULT VALUE OF YOUR PARAMETERS).
	// ADD MORE LINES IF NECESSARY (YOU SHOULD HAVE nbparamstat[0] LINES).
	paramstat[0] = 0.0;
     }
// THE FOLLOWING 8 LINES SHOULD NOT BE MODIFIED.
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i=j;i<50;i++) name[i][0] = space[0];
      return;
    }
	
// MODIFY THESE LINES TO SUIT YOUR NEEDS. SEE OTHER statxxx.cpp FOR MORE EXAMPLES.
    double mu;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      mu = 0.0;
      paramstat[0] = 0.0;
    } else if (nbparamstat[0] == 1) {
      mu = paramstat[0];
    } else {
      return;
    }


    if (n>3) {
// Computation of the value of the test statistic
      double stat;
//    YOU NEED TO COMPUTE THE VALUE OF YOUR TEST STATISTIC (VARIABLE stat) USING 
//    THE x[i] (i=0,..,n-1) AND THE VALUES OF THE PARAMETERS (E.G., mu AS DEFINED ABOVE).
//  ADD YOUR OWN CODE HERE.	
// THE LAST LINE SHOULD BE:
    statistic[0] = stat;

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
   #include "pvalues/pvaluej.cpp"  // YOU HAVE TO MODIFY THE CONTENTS OF THE FILE pvaluej.cpp TO SUIT YOUR NEEDS.
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      // YOU MAY NEED TO MODIFY THESE LINES DEPENDING ON THE VALUES OF alter IN YOUR CASE. SEE OTHER statxxx.cpp FOR MORE EXAMPLES.
	  if (usecrit[0] == 1) { // We use the provided critical values
	    if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	    } else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	    } else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
	  } else {
	    if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
	  }
    }
	
    
// IF APPLICABLE, FREE UNUSED ARRAYS OF POINTERS YOU MAY HAVE DEFINED ABOVE.

}

// We return
    return;
   
        
  }
  
  
  
}
