// Title: Nout(a)
// Ref. (book or article): Romão, Xavier, Delgado, Raimundo and Costa, Aníbal, 
//			   "An empirical power comparison of univariate goodness-of-fit tests for normality",
//	                   J. Stat. Comput. Simul. 80 (2010), no. 5-6, 545--591.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law33(int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Nout(a)$";
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
    double a;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      a = 1.0;
      params[0] = 1.0;
    } else if (nbparams[0] == 1) {
      a = params[0];
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (a != 1.0 && a != 2.0 && a != 3.0 && a != 4.0 && a != 5.0) {
      warning("a should be a value in {1,2,3,4,5} in law33!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double rf_rnorm(double mu, double sigma);
    double xq1, xq3, IQR;
    xq1 = -0.6744897502;
    xq3 = 0.6744897502;
    IQR = 1.348979500;
    if (((int)a != 1) && ((int)a != 2) && ((int)a != 3) && ((int)a != 4) && ((int)a != 5)) return; // protection
    if ((int)a == 1) {
      x[0] = xq3 + 2.0 * IQR;
      for (i = 1; i < n; i++) x[i] = Rf_rnorm(0.0, 1.0);
    }
    if ((int)a == 2) {
      x[0] = xq3 + 3.0 * IQR;
      for (i = 1; i < n; i++) x[i] = Rf_rnorm(0.0, 1.0);
    }
    if ((int)a == 3) {
      x[0] = xq3 + 2.0 * IQR;
      x[1] = xq3 + 3.0 * IQR;
      for (i = 2; i < n; i++) x[i] = Rf_rnorm(0.0, 1.0);
    }
    if ((int)a == 4) {
      x[0] = xq3 + 2.0 * IQR;
      x[1] = xq1 - 2.0 * IQR;
      for (i = 2; i < n; i++) x[i] = Rf_rnorm(0.0, 1.0);
    }
    if ((int)a == 5) {
      x[0] = xq3 + 2.0 * IQR;
      x[1] = xq1 - 2.0 * IQR;
      x[2] = xq3 + 3.0 * IQR;
      x[3] = xq1 - 3.0 * IQR;
      for (i = 4; i < n; i++) x[i] = Rf_rnorm(0.0, 1.0);
    }
    if (setseed[0] == 1) PutRNGstate();
    
// If applicable, we free the unused array of pointers. Then we return.
    return;
    
  }
  
}
