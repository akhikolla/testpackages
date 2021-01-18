#include <Rmath.h>

void kern1Gauss(int *n,             /* number of u               */
		double *u,          /* distances (positive)      */
		double *range,      /* range parameter           */
		double *var,        /* variance                  */
		double *kern,       /* n kernel values           */ 
		double *dkern) {    /* n kernel gradient values  */  
  double z, z2;

  for (int i = 0; i < *n; i++) {
    z = u[i] / range[0];
    z2 = z * z;
    kern[i] = exp(-z2 / 2.0) * *var;
    dkern[i] = -z2 * kern[i] / range[0]; 
  }
}

void kern1Exp(int *n,             /* number of u               */
	      double *u,          /* distances (positive)      */
	      double *range,      /* range parameter           */
	      double *var,        /* variance                  */
	      double *kern,       /* n kernel values           */ 
	      double *dkern) {    /* n kernel gradient values  */
  double z;
  
  for (int i = 0; i < *n; i++) {
    z = u[i] / range[0];
    kern[i] = exp(-z) * *var;
    dkern[i] = -z * kern[i] / range[0]; 
  }
}

void kern1PowExp(int *n,                /* number of u               */
		 double *u,             /* distances (positive)      */
		 double *rangeShape,    /* range & shape vector      */
		 double *var,           /* variance                  */
		 double *kern,          /* n kernel values           */ 
		 double *dkern) {       /* n kernel gradient values  */
  double z, zAlpha;
		 
  for (int i = 0; i < *n; i++) {
    z = u[i] / rangeShape[0];
    zAlpha = pow(z, rangeShape[1]); 
    kern[i] = exp(-zAlpha) * *var;
    dkern[i] =  z * kern[i] * rangeShape[1] / rangeShape[0];
    if (z < 1e-6) {
      dkern[i + *n] = 0.0;
    } 
    else {
      dkern[i + *n] = -log(z) * zAlpha * kern[i];
    }
  }
}
