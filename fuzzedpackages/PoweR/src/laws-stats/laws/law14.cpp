// Title: Mean of k r.v. following a Uniform distribution 
// Ref. (book or article): Quesenberry and Miller, "Power studies of some tests for uniformity", 
//                         Journal of Statistical Computation and Simulation (1977), 5:3, 169--191 see p. 179
//	                   Bates, Grace E., "Joint distributions of time intervals for the occurrence of successive accidents in a generalized Pólya scheme.",
//	                   Ann. Math. Statist. 26 (1955), 705--720. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law14 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$AveUnif(k,a,b)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 2;
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
    double k;
    double a, b;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      k = 2;
      a = 0.0;
      b = 1.0;
      params[0] = 2;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      k = params[0];
      a = 0.0;
      b = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      k = params[0];
      a = params[1];
      b = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 3) {
      k = params[0];
      a = params[1];
      b = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }

// If necessary, we check if some parameter values are out of parameter space
    double intpart;
    if (k < 1.0 || modf(k, &intpart) != 0.0 || a > b) {
      warning("k should be a positive integer value and a should not be > b in law14!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    for (i = 0; i < n; i++) {
      x[i] = Rf_runif(a,b);
      for (j = 1; j <= ((int)k - 1); j++) {
	x[i] = x[i] + Rf_runif(a, b);
      }
      x[i] = x[i] / k;
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
