// Title: S(a,b)
// Ref. (book or article): Chambers, J. M., Mallows, C. L. and Stuck, B. W. (1976),
//						   "A method for simulating stable random variables", 
//						   J. Amer. Statist. Assoc. 71 (1976), no. 354, 340--344.  
//						   or see function rstable() in package stabledist

#include <R.h>
#include "Rmath.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

extern "C" {

  void law25 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Stable(\\alpha,\\beta,c,\\mu)$";
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
    double alpha, beta, c, mu;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      alpha = 2.0;
      beta = 0.0;
      c = 1.0;
      mu = 0.0;
      params[0] = 1.0;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = 0.0;
      c = 1.0;
      mu = 0.0;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = params[1];
      c = 1.0;
      mu = 0.0;
      params[2] = 1.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      alpha = params[0];
      beta = params[1];
      c = params[2];
      mu = 0.0;
      params[3] = 0.0;
    } else if  (nbparams[0] == 4) {
      alpha = params[0];
      beta = params[1];
      c = params[2];
      mu = params[3];
    } else {
      error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (alpha <= 0.0 || alpha > 2.0 || beta < -1.0 || beta > 1.0 || c <= 0.0) {
      warning("Some parameter(s) value(s) are out of parameter space in in law25!\nCorrect values are: alpha>0.0 && alpha<=2.0 && beta>=-1.0 && beta<=1.0 && c>0.0\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();	
    double Rf_rcauchy(double location, double scale);
    double Rf_runif(double a, double b);
    double *U;
    U = new double [n];
    double pi2, btanpa, theta0, tmp = 0.0, cc, theta = 0.0, w = 0.0, atht;
    pi2 = M_PI / 2.0;
    if (alpha == 1.0 && beta == 0.0) {
      for (i = 0; i < n; i++) {
	U[i] = Rf_rcauchy(0.0, 1.0);
      }
        // Otherwise, if alpha is different from 1:
    } else {
      btanpa = beta * tan(pi2 * alpha);
      tmp    = MAX(-pi2, atan(btanpa) / alpha);
      theta0 = MIN(tmp, pi2);
      cc     = R_pow(1.0 + R_pow(btanpa, 2.0), 1.0 / (2.0 * alpha));		
      for (i = 0; i < n; i++) {
	theta = M_PI * (Rf_runif(0.0, 1.0) - 1.0 / 2.0);
	w = -log(Rf_runif(0.0, 1.0));
	atht = alpha * (theta + theta0);
	U[i] = (cc * sin(atht) / R_pow(cos(theta), 1.0 / alpha)) * R_pow(cos(theta - atht) / w, (1.0 - alpha) / alpha) - btanpa;
      }
    }	
    for (i = 0; i < n; i++) x[i] = U[i] * c + mu;
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;

// We return
    return;
 
  }
  
}



