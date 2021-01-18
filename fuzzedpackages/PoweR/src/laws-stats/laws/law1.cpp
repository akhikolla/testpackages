// Title: Laplace(mu,b)
// Ref. (book or article): See function urlaplace() from Runuran package

#include <R.h>
#include "Rmath.h"

extern "C" {
  
  void law1 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Laplace(\\mu,b)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
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
    double mu, b;
    if (nbparams[0] == 0) {
      nbparams[0] = 2;
      mu = 0.0;
      b = 1.0;
      params[0] = 0.0;
      params[1] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 2;
      mu = params[0];
      b = 1.0;
      params[1] = 1.0;
    } else if  (nbparams[0] == 2) {
      mu = params[0];
      b = params[1];
    } else {
      error("Number of parameters should be at most: 2");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (b <= 0.0) {
      warning("b should not be <= 0 in law1!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);  
    double *U;
    U = new double [n];
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0, 1.0);
    for (i = 0; i < n; i++)   {
      if (U[i] > 0.5) x[i] = mu - b * log(2.0 * (1.0 - U[i])); 
	  else x[i] = mu + b * log(2.0 * U[i]);
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
