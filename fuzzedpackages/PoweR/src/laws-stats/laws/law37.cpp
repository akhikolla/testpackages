// Title: Normal-inverse Gaussian(alpha,beta,mu,delta)
// Ref. (book or article): Atkinson, A. C. (1982), "The simulation of generalized inverse Gaussian and hyperbolic random variables",
//	                       SIAM J. Sci. Statist. Comput. 3 (1982), No. 4, 502--515. 
//                         or see function rnig() in package fBasics

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law37 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
   const char *nom = "$NIG(\\alpha,\\beta,\\delta,\\mu)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 4;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
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
    double alpha, beta, mu, delta;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      alpha = 1.0;
      beta = 0.0;
      delta = 1.0;
      mu = 0.0;
      params[0] = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = 0.0;
      delta = 1.0;
      mu = 0.0;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = params[1];
      delta = 1.0;
      mu = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = params[1];
      delta = params[2];
      mu = 0.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 4) {
      alpha = params[0];
      beta = params[1];
      delta = params[2];
      mu = params[3];
    } else {
      error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (alpha <= 0.0 || delta < 0.0 || beta < -alpha || beta > alpha) {
      warning("correct values of the parameters are: alpha > 0 && delta >= 0 && beta in (-alpha,alpha) in law37!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate(); 	
    double Rf_runif(double a, double b);
    double Rf_rnorm(double a, double b);
    double z1(double v, double d, double g);
    double z2(double v, double d, double g);
    double pz1(double v, double d, double g);
    double *U;
    U = new double [n];    
    double *V;
    V = new double [n];
    double *Z;
    Z = new double [n];
    double *S;
    S = new double [n];
    double gamma;
    // Settings:
    gamma = sqrt(alpha * alpha - beta * beta);
    // GAMMA:
    if (gamma == 0) {
        // GAMMA = 0:
      for (i=0;i<n;i++) {
	V[i] = R_pow(rnorm(0.0,1.0),2.0);
	Z[i] = delta*delta/V[i];
	x[i] = sqrt(Z[i])*rnorm(0.0,1.0);
      }
    } else {
        // GAMMA > 0:
      for (i = 0; i < n; i++) {
	U[i] = Rf_runif(0.0, 1.0);
	V[i] = R_pow(Rf_rnorm(0.0, 1.0), 2.0);
	S[i] = U[i] - pz1(V[i], delta, gamma);
	if (S[i] < 0) {
	  S[i] = 1.0;
	} else {
	  if (S[i] == 0) {
	    S[i] = 0.5;
	  } else {
	    if (S[i] > 0) {
	      S[i] = 0.0;
	    }
	  }
	}
	Z[i] = z1(V[i], delta, gamma) * S[i] + z2(V[i], delta, gamma) * (1.0 - S[i]);
	x[i] = mu + beta * Z[i] + sqrt(Z[i]) * Rf_rnorm(0.0, 1.0);
      }
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    delete[] V;
    delete[] Z;
    delete[] S;	
	
// We return
    return;
    
  }
	
	
  // function z1	
  double z1(double v, double d, double g) {
	return d/g + v/(2.0*g*g) - sqrt(v*d/(g*g*g) + v*v/(4.0*g*g*g*g));
  }

  // function z2
  double z2(double v, double d, double g) { 
	return (d*d/(g*g))/z1(v,d,g);
  }
       
  // function pz1
  double pz1(double v, double d, double g) { 
	return d/(d + g * z1(v,d,g));
  }
  
  
}
