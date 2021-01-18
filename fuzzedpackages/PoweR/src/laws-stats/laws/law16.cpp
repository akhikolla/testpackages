// Title: Laws in V: (Y+0.5) if Y<0.5 and (Y-0.5) if Y>0.5 with Y a mean of (j+1) independent Uniform r.v.
// Ref. (book or article): Quesenberry and Miller, "Power studies of some tests for uniformity", 
//                         Journal of Statistical Computation and Simulation (1977), 5:3, 169--191 see p. 179

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law16 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$VUnif(j)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 1.0;
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
	// nbunif = params[0] + 1 ALWAYS !!!
    double nbunif;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      nbunif = 2.0;
      params[0] = 1.0;
    } else if (nbparams[0] == 1) {
      nbunif = params[0] + 1.0;
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    double intpart;
    if (nbunif <= 1.0 || modf(nbunif, &intpart) != 0.0) {
      warning("nbunif should be a positive integer value in law16!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double temp;
    int k;
    temp = 0.0;
    for (i = 0; i < n; i++) {
      temp = Rf_runif(0.0, 1.0);
      for (k = 1; k < (int)nbunif; k++) {
	temp = temp + Rf_runif(0.0, 1.0);
      }
      temp = temp / (double)nbunif;
      if (temp < 0.5) 
	x[i] = temp + 0.5 ;
      else 
	x[i] = temp - 0.5 ;
    }
    if (setseed[0] == 1) PutRNGstate();	
    
// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
