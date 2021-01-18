// Title: Skew Normal(xi,omega,alpha)
// Ref. (book or article): Azzalini, Adelchi(2005), "The skew-normal distribution and related multivariate families",
//						   With discussion by Marc G. Genton and a rejoinder by the author, 
//                         Scand. J. Statist. 32 (2005), no. 2, 159--200. 
//                         or see function rsnorm() in package VGAM

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law21 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$SkewN(\\xi,\\omega,\\alpha)$";
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
    double xi, omega, alpha;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      xi = 0.0;
      omega = 1.0;
      alpha = 0.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      xi = params[0];
      omega = 1.0;
      alpha = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      xi = params[0];
      omega = params[1];
      alpha = 0.0;
      params[2] = 0.0;
    } else if  (nbparams[0] == 3) {
      xi = params[0];
      omega = params[1];
      alpha = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }
   
// If necessary, we check if some parameter values are out of parameter space
    if (omega <= 0.0) {
      warning("omega should not be <= 0 in law21!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rnorm(double mean, double sd);
    double delta;
    delta = alpha / sqrt(1.0 + alpha * alpha);
    double *U0;
    U0 = new double [n];
    for (i = 0; i < n; i++) U0[i] = Rf_rnorm(0.0,1.0);
    for (i = 0; i < n; i++)   {
      if (U0[i]>=0) x[i] = xi + omega * (delta * U0[i] + sqrt(1.0 - delta * delta) * Rf_rnorm(0.0, 1.0)); 
      else x[i] = xi + omega * (-(delta * U0[i] + sqrt(1.0 - delta * delta) * Rf_rnorm(0.0, 1.0))); 
    }
    if (setseed[0] == 1) PutRNGstate();
    
// If applicable, we free the unused array of pointers. Then we return.
    delete[] U0;
    return;
    
  }
  
}
