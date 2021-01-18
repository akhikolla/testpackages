// Title: SU^{j+1}+(1-S)(1-U^{j+1}) with Pr(S=0)=Pr(S=1)=(1/2)
// Ref. (book or article): Quesenberry and Miller, "Power studies of some tests for uniformity", 
//                         Journal of Statistical Computation and Simulation (1977), 5:3, 169--191 see p. 179

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law15 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$UUnif(j)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 1.0;
     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i = j; i < 50; i++) name[i][0] = space[0];
      return;
    }

// Initialization of the parameters
    double power;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      power = 2.0;
      params[0] = 1.0;
    } else if (nbparams[0] == 1) {
      power = params[0] + 1.0 ;
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (power < 0.0) {
      warning("power should not be < 0 in law15!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate(); 
    double val;  
    double Rf_runif(double a, double b);
    for (i = 0; i < n; i++) {
      val = Rf_runif(0.0, 1.0);
      if (Rf_runif(0.0, 1.0) < 0.5) 
	x[i] = R_pow(val, power);
      else 
	x[i] = 1.0 - R_pow(val, power);
    }
    if (setseed[0] == 1) PutRNGstate();
    
// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
