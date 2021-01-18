// Title: Generalized Pareto GP(mu,sigma,xi)
// Ref. (book or article): Choulakian, V.; Stephens, M. A. Goodness-of-fit tests for the generalized Pareto distribution.  Technometrics  43  (2001),  no. 4, 478–484.
//                         or see function rgpd() in package POT

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law23 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$GP(\\mu,\\sigma,\\xi)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
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
    double mu, sigma, xi;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      mu = 0.0;
      sigma = 1.0;
      xi = 0.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = 1.0;
      xi = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = params[1];
      xi = 0.0;
      params[2] = 0.0;
    } else if  (nbparams[0] == 3) {
      mu = params[0];
      sigma = params[1];
      xi = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }
   
// If necessary, we check if some parameter values are out of parameter space
    if (sigma <= 0.0) {
      warning("sigma should not be <= 0 in law23!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rexp(double lambda);
    double *U;
    U = new double [n];	
    if (xi == 0.0) {
      for (i = 0; i < n; i++) {
        U[i] = Rf_rexp(1.0);
        x[i] = mu + sigma * U[i]; 
      }	  
    }
    else {
      for (i = 0; i < n; i++) {
        U[i] = Rf_runif(0.0, 1.0);
        x[i] = mu + sigma * (1.0 - R_pow(U[i], xi)) / xi; 
      }
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
