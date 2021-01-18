#include <math.h>

double C_covWhiteNoise(const double *x1, const int *n1, 
                       const double *x2, const int *n2, 
                       const int *d, 
		       const int *i1, const int *i2, 
                       const double *var) {
  double s = 0.;
  for (int k = 0; k < *d; k++) {
    s += fabs((x1[*i1 + *n1 * k] - x2[*i2 + *n2 * k]));
  }
  if (s < 0.000000000000001) return(*var);
  else return(0.);
}

void C_covMat1Mat2_WhiteNoise(const double *x1, const int *n1, 
                              const double *x2, const int *n2, 
                              const int *d, const double *var, 
                              double *ans) {
    
  for (int i1 = 0; i1 < *n1; i1++) {
    for (int i2 = 0; i2 < *n2; i2++) {				
      ans[i1 + *n1 * i2] = C_covWhiteNoise(x1, n1, x2, n2, d, &i1, &i2, var);
    }
  }
  
} 
