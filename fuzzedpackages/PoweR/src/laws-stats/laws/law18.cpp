// Title: Symmetrical Tukey(l)
// Ref. (book or article):  Joiner, Brian L.; Rosenblatt, Joan R. (1971), "Some Properties of the Range in Samples from Tukey's Symmetric Lambda Distributions", 
//                          Journal of the American Statistical Association 66 (334): 394–399, 


#include <R.h>
#include "Rmath.h"

extern "C" {

  void law18 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Tukey(l)$";
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
    double l;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      l = 1.0;
      params[0] = 1.0;
    } else if (nbparams[0] == 1) {
      l = params[0];
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_runif(double a, double b);
    double *U;
    U = new double [n];
    for (i = 0; i < n; i++) U[i] = Rf_runif(0.0,1.0);
    if (l == 0) {
       for (i = 0; i < n;i++)   {
	 x[i] = log(U[i]) - log(1.0 - U[i]);
      }     
    } else {
      for (i = 0; i < n; i++)   {
	x[i] = (R_pow(U[i], l) - R_pow(1 - U[i], l)) / l;
      }
    }
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] U;
    return;
    
  }
  
}
