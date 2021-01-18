// Title: LoConN(p,m)
// Ref. (book or article): Kallenberg, Wilbert C. M. and Ledwina, Teresa (1997), 
//						   "Data driven smooth tests for composite hypotheses: comparison of powers", 
//						   J. Statist. Comput. Simulation 59 (1997), No. 2, 101--121. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law19 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$LoConN(p,m)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.2;
      params[1] = 3.0;
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
    double p, m;
    if (nbparams[0] == 0) {
      nbparams[0] = 2;
      p = 0.2;
      m = 3.0;
      params[0] = 0.2;
      params[1] = 3.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 2;
      p = params[0];
      m = 3.0;
      params[1] = 0.2;
    } else if  (nbparams[0] == 2) {
      p = params[0];
      m = params[1];
    } else {
      error("Number of parameters should be at most: 2");
    }
   
// If necessary, we check if some parameter values are out of parameter space
    if (p < 0.0 || p > 1.0) {
      warning("p should be in [0,1] in law19!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double Rf_rnorm(double mean, double sd);
    double *U;
    U = new double [n];
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0,1.0);
    for (i = 0; i < n; i++)   {
      if (U[i] < p) x[i] = Rf_rnorm(m, 1.0); else x[i] = Rf_rnorm(0.0, 1.0);
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
