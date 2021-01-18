// Title: Log-Pareto-tail-normal(alpha, mu, sigma, log.d)
// Ref. (book or article): Desgagné, Alain. Robustness to outliers in location–scale parameter model using log-regularly varying distributions. Ann. Statist. 43 (2015), no. 4, 1568--1595. doi:10.1214/15-AOS1316. http://projecteuclid.org/euclid.aos/1434546215.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law40 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$LPtn(alpha,mu,sigma)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 1.959964;
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
    double alpha, mu, sigma;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      alpha = 1.959964;
      mu = 0.0;
      sigma = 1.0;
      params[0] = 1.959964;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      alpha = params[0];
      mu = 0.0;
      sigma = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      alpha = params[0];
      mu = params[1];
      sigma = 1.0;
      params[2] = 1.0;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 3;
      alpha = params[0];
      mu = params[1];
      sigma = params[2];
   } else {
      error("Number of parameters should be at most: 3");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (sigma <= 0.0 || alpha <= 1.0) {
      warning("sigma should be > 0 and alpha should be larger than 1 in law40!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rbinom(double nin, double pp);
    double Rf_qnorm5(double p, double mean, double sd, int lower_tail, int log_p);
    double Rf_pnorm5(double q, double mean, double sd, int lower_tail, int log_p);
    double Rf_dnorm4(double x, double mean, double sd, int log);  
    double q, beta, n1, n2, n3, Galpha;
    double *U, *U1, *U2, *U3, *x1, *x2, *x3;
    U = new double [n];
    q = 1.0 - 2.0 * Rf_pnorm5(alpha, 0.0, 1.0, 0, 0);
    beta = 1.0 + (2.0 * Rf_dnorm4(alpha, 0.0, 1.0, 0) * alpha * log(alpha)) / (1.0 - q);
    n1 = Rf_rbinom((double)n, q);
    n2 = Rf_rbinom((double)n - n1, 0.5);
    n3 = (double)n - n1 - n2;  
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0, 1.0);
    U1 = new double [(int)n1];
    for (i = 0; i < (int)n1; i++) U1[i] = U[i]; 

    U2 = new double [(int)n2];
    for (i = 0; i < (int)n2; i++) U2[i] = U[(int)n1 + i]; 
    U3 = new double [(int)n3];
    for (i = 0; i < (int)n3; i++) U3[i] = U[(int)n1 + (int)n2 + i]; 
    Galpha = (1.0 + q) / 2.0;
    x1 = new double [(int)n1];
    for (i = 0; i < (int)n1; i++) x1[i] = Rf_qnorm5((2.0 * Galpha - 1.0) * U1[i] + 1.0 - Galpha, 0.0, 1.0, 1, 0);
    x2 = new double [(int)n2];
    for (i = 0; i < (int)n2; i++) x2[i] = exp(log(alpha) / R_pow(U2[i], 1.0 / (beta - 1.0)));
    x3 = new double [(int)n3];
    for (i = 0; i < (int)n3; i++) x3[i] = -exp(log(alpha) / R_pow(U3[i], 1.0 / (beta - 1.0)));
    for (i = 0; i < (int)n1; i++) x[i] = mu + sigma * x1[i];
    for (i = 0; i < (int)n2; i++) x[i + (int)n1] = mu + sigma * x2[i];
    for (i = 0; i < (int)n3; i++) x[i + ((int)n1 + (int)n2)] = mu + sigma * x3[i];
    
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    delete[] U1;
    delete[] x1;
    delete[] U2;
    delete[] x2;
    delete[] U3;
    delete[] x3;

  return;
    
  }
  
}
