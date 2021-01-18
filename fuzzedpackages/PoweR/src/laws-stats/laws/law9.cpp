// Title: Chi2(df)
// Ref. (book or article): See function rchisq() from package stats

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law9 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Chi-squared(df)$";
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
    double df;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      df = 1.0;
      params[0] = 1.0;
    } else if (nbparams[0] == 1) {
      df = params[0];
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (df < 0.0) {
      warning("df should not be < 0 in law9!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_rchisq(double df);
    for (i = 0; i < n; i++) x[i] = Rf_rchisq(df);
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
  }
  
}

