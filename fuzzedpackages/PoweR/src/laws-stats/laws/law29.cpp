// Title: Generalized Arcsine(alpha)
// Ref. (book or article): Feller, William (1968), "An introduction to probability theory and its applications. Vol. I.", 
//                         Third edition John Wiley & Sons, Inc., New York-London-Sydney 1968 xviii+509 pp.
//                         Feller, William (1971), "An introduction to probability theory and its applications. Vol. II.", 
//                         Second edition John Wiley & Sons, Inc., New York-London-Sydney 1971 xxiv+669 pp.   
//                         MAXIMILIAN THALER, "THE DYNKIN-LAMPERTI ARC-SINE LAWS FOR MEASURE PRESERVING TRANSFORMATIONS", TRANSACTIONS OF THE AMERICAN MATHEMATICAL SOCIETY
//                         Volume 350, Number 11, November 1998, Pages 4593–4607 S 0002-9947(98)02312-5 See page 4594 (with \alpha --> -\alpha)

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law29 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$GArcSine(\\alpha)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.5;
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
    double alpha;
    if (nbparams[0] == 0) {
      nbparams[0] = 1;
      alpha = 0.5;
      params[0] = 0.5;
    } else if (nbparams[0] == 1) {
      alpha = params[0];
    } else {
      error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (alpha <= 0.0 || alpha >= 1.0) {
      warning("You should take 0 < alpha <1 in law29!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    double Rf_rbeta(double a, double b);
    for (i = 0; i < n; i++) x[i] = Rf_rbeta(1.0 - alpha, alpha);
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
 }

}



