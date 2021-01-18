// Title: Frechet(mu,sigma,alpha)
// Ref. (book or article): Castillo, Enrique, H., Ali S. and Balakrishnan, N. and Sarabia, J. M. (2005),
//                         "Extreme value and related models with applications in engineering and science",
//                         Wiley Series in Probability and Statistics. Wiley-Interscience [John Wiley & Sons], Hoboken, NJ, 2005. xiv+362 pp. ISBN: 0-471-67172-X page 198
//                         or see function rfrechet() in package VGAM  

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law27 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Frechet(\\mu,\\sigma,\\alpha)$";
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
    double mu, sigma, alpha;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      mu = 0.0;
      sigma = 1.0;
      alpha = 1.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = 1.0;
      alpha = 1.0;
      params[1] = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      mu = params[0];
      sigma = params[1];
      alpha = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 3) {
      mu = params[0];
      sigma = params[1];
      alpha = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (sigma <= 0.0 || alpha <= 0.0) {
      warning("sigma and alpha should be > 0 in law27!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double Rf_runif(double a, double b);
    for (i = 0; i < n; i++) x[i] = mu + sigma * R_pow(-log(Rf_runif(0.0, 1.0)), -1.0 / alpha);
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
 }

}



