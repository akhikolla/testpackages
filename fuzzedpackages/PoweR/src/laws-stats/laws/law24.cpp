// Title: GED(mu,sigma,p)
// Ref. (book or article): Nadarajah, Saralees (2005), "A generalized normal distribution",
//						   J. Appl. Stat. 32 (2005), no. 7, 685--694.  
//                         See also function rnormp() in package normalp but note the reparameterization: \sigma_p --> \sigma/(p^{1/p})

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law24 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$GED(\\mu,\\sigma,p)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
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
    double mu, sigma, p;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      mu = 0.0;
      sigma = 1.0;
      p = 1.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = 1.0;
      p = 1.0;
      params[1] = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = params[1];
      p = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 3) {
      mu = params[0];
      sigma = params[1];
      p = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (sigma <= 0.0 || p <= 0.0) {
      warning("p or sigma should not be <= 0 in law24!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double Rf_runif(double a, double b);
    double shape, scale;
    shape = 1 / p;
    scale = p;
    double *U;
    U = new double [n];
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0, 1.0);
    double Rf_rgamma(double shape, double scale);
    for (i = 0; i < n; i++) {
      if (U[i] > 0.5) x[i] = mu + sigma * R_pow(Rf_rgamma(shape, scale) / p, 1.0 / p); else x[i] = mu - sigma * R_pow(Rf_rgamma(shape, scale) / p, 1.0 / p);
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
