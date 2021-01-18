// Title: modified Asymmetric Power Distribution APD(theta,phi,alpha,lambda)
// Ref. (book or article): Our second paper with Desgagne and Leblanc. And see also Komunjer, Ivana "Asymmetric Power Distribution: Theory and Applications to Risk Measurement",  
//                         Journal of Applied Econometrics (2007), 22, 891-921.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law39 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$modAPD(mu,sigma,theta1,theta2)$";
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
    double theta1, theta2, mu, sigma;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      mu = 0.0;
      sigma = 1.0;
      theta1 = 0.5;
      theta2 = 2.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      mu = params[0];
      sigma = 1.0;
      theta1 = 0.5;
      theta2 = 2.0;
      params[1] = 1.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      mu = params[0];
      sigma = params[1];
      theta1 = 0.5;
      theta2 = 2.0;
      params[2] = 0.5;
      params[3] = 2.0;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      mu = params[0];
      sigma = params[1];
      theta1 = params[2];
      theta2 = 2.0;
      params[3] = 2.0;
    } else if  (nbparams[0] == 4) {
      mu = params[0];
      sigma = params[1];
      theta1 = params[2];
      theta2 = params[3];
   } else {
      error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (sigma <= 0.0 || theta1 <= 0.0 || theta1 >= 1.0 || theta2 <= 0.0) {
      warning("theta2 and sigma should be > 0 and you must take 0 < theta1 < 1 in law39!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rgamma(double shape, double scale);
    double invtheta2 = 1.0 / theta2;
    double twopowinvtheta2 = R_pow(2.0, invtheta2);
    double delta = (2.0 * R_pow(theta1, theta2) * R_pow(1.0 - theta1, theta2)) / (R_pow(theta1, theta2) + R_pow(1.0 - theta1, theta2));
    for (i = 0; i < n; i++) {
      if (Rf_runif(0.0, 1.0) <= theta1) {
	x[i] = mu + sigma * twopowinvtheta2 * (-theta1 * R_pow(Rf_rgamma(invtheta2, 1.0) / delta, invtheta2)); // Petit doute si le 2 devrait etre au debut: 2 * (theta + phi ...)
      } else {
	x[i] = mu + sigma * twopowinvtheta2 * ((1.0 - theta1) * R_pow(Rf_rgamma(invtheta2, 1.0) / delta, invtheta2)); // Petit doute si le 2 devrait etre au debut: 2 * (theta + phi ...)
      }
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
