// Title: MixN(p,m)
// Ref. (book or article): Romão, Xavier, Delgado, Raimundo and Costa, Aníbal, 
//			   "An empirical power comparison of univariate goodness-of-fit tests for normality",
//                         J. Stat. Comput. Simul. 80 (2010), no. 5-6, 545--591.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law31(int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$MixN(p,m,d)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.5;
      params[1] = 0.0;
      params[2] = 1.0;
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
    double p, m, d;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      p = 0.5;
      m = 0.0;
      d = 1.0;
      params[0] = 0.5;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      p = params[0];
      m = 0.0;
      d = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      p = params[0];
      m = params[1];
      d = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 3) {
      p = params[0];
      m = params[1];
      d = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }
   
// If necessary, we check if some parameter values are out of parameter space
    if (d <= 0.0 || p < 0.0 || p > 1.0) {
      warning("d should be >0 and p should be in [0,1] in law31!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rnorm(double mean, double sd);
    double *U;
    U = new double [n];
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0,1.0);
    for (i = 0; i < n;i++)   {
      if (U[i] < p) x[i] = Rf_rnorm(m, d); else x[i] = Rf_rnorm(0.0, 1.0);
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
