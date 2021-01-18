// Title: Asymmetric Power Distribution APD(theta,phi,alpha,lambda)
// Ref. (book or article): Komunjer, Ivana "Asymmetric Power Distribution: Theory and Applications to Risk Measurement",  
//                         Journal of Applied Econometrics (2007), 22, 891-921.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law38 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$APD(theta,phi,alpha,lambda)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 4;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.5;
      params[3] = 2.0;
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
    double alpha, lambda, theta, phi;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      theta = 0.0;
      phi = 1.0;
      alpha = 0.5;
      lambda = 2.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      theta = params[0];
      phi = 1.0;
      alpha = 0.5;
      lambda = 2.0;
      params[1] = 1.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      theta = params[0];
      phi = params[1];
      alpha = 0.5;
      lambda = 2.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      theta = params[0];
      phi = params[1];
      alpha = params[2];
      lambda = 2.0;
      params[3] = 2.0;
    } else if  (nbparams[0] == 4) {
      theta = params[0];
      phi = params[1];
      alpha = params[2];
      lambda = params[3];
   } else {
      error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (phi <= 0.0 || alpha <= 0.0 || alpha >= 1.0 || lambda <= 0.0) {
      warning("lambda and phi should be > 0 and you must take 0 < alpha < 1 in law38!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rgamma(double shape, double scale);
    double delta = (2.0 * R_pow(alpha, lambda) * R_pow(1.0 - alpha, lambda)) / (R_pow(alpha, lambda) + R_pow(1.0 - alpha, lambda));
    for (i = 0; i < n; i++) {
      if (Rf_runif(0.0, 1.0) <= alpha) {
	x[i] = theta + phi * (-alpha * R_pow(Rf_rgamma(1.0 / lambda, 1.0) / delta, 1.0 / lambda));
      } else {
	x[i] = theta + phi * ((1.0 - alpha) * R_pow(Rf_rgamma(1.0 / lambda, 1.0) / delta, 1.0 / lambda));
      }
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
