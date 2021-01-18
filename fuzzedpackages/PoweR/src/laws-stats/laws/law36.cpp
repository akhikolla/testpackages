// Title: Asymmetric Laplace(mu,b,k)
// Ref. (book or article): Kotz, Samuel, Kozubowski, Tomasz J. and Podgórski, Krzysztof (2001),
// 	                       "The Laplace distribution and generalizations. A revisit with applications to communications, economics, engineering, and finance." 
//						   Birkhäuser Boston, Inc., Boston, MA, 2001. xviii+349 pp. ISBN: 0-8176-4166-1, see Equation (3.0.8) page 134
//			   See also function ralap() from package VGAM.
//                         Note that this is a reparameterization: \sigma --> \sigma/\sqrt{2}

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law36 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
    const char *nom = "$ALaplace(\\mu,b,k)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 2.0;
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
    double mu, b, k;
    if (nbparams[0] == 0) {
      nbparams[0] = 3;
      mu = 0.0;
      b = 1.0;
      k = 2.0;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 2.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 3;
      mu = params[0];
      b = 1.0;
      k = 2.0;
      params[1] = 1.0;
      params[2] = 2.0;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 3;
      mu = params[0];
      b = params[1];
      k = 2.0;
      params[2] = 2.0;
    } else if  (nbparams[0] == 3) {
      mu = params[0];
      b = params[1];
      k = params[2];
    } else {
      error("Number of parameters should be at most: 3");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (b < 0.0 || k <= 0.0) {
      warning("b should be >=0 and k should be >0 in law36!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);  
    for (i = 0; i < n; i++) {
      x[i] = mu + b * log(R_pow(Rf_runif(0.0, 1.0), k) / R_pow(Rf_runif(0.0, 1.0), 1.0 / k)) / sqrt(2.0);
    }
    if (setseed[0] == 1) PutRNGstate();		
    
// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
