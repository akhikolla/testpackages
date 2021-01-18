// Title: Statistique de test de Epps-Pulley with our new volcano weight involving the location parameter mu
// Ref. (book or article): Epps, T.W. and Pulley, L.B. (1983), A test of normality based on empirical characteristic function, 
//						   Biometrika, Vol. 70, No. 3, pp. 723-726.

#include <R.h>
#include "Rmath.h"

extern "C" {
  //#include "libraries/Faddeeva.cc"

  void stat90(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T^{LS}_4(\\alpha,\\mu)$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 1.15;
	paramstat[1] = 1.15;
     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i=j;i<50;i++) name[i][0] = space[0];
      return;
    }
	
// Initialization of the parameters
    double mu, alpha;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      mu = 1.15;
      alpha = 1.15;
      paramstat[0] = 1.15;
      paramstat[1] = 1.15;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      mu = paramstat[0];
      alpha = 1.15;
      paramstat[1] = 1.15;
    } else if (nbparamstat[0] == 2) {
      mu = paramstat[0];
      alpha = paramstat[1];
    } else {
      return;
    }

	
    if (n>3) {
// Computation of the value of the test statistic
//      double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
      // //      extern double complex Faddeeva_erf(double complex z, double relerr);
      double statTEP, term1 = 0.0, term2 = 0.0, term3, xbar = 0.0, S2 = 0.0, S, tmp1, tmp2, alpha2, mu2, xx, xx2, mux, cosmux, sinmux;
      double y, y2, sinmuy, cosmuy, calphamu;
      int k;
      cmplx z;

	  // calculate xbar
      for (i = 0; i < n; i++) {
	xbar = xbar + x[i];
      }
      xbar = xbar / (double)n;
      
	  // calculate S^2
      for (i = 0; i < n; i++) {
	S2 = S2 + R_pow(x[i] - xbar, 2.0);
      }
      S2 = S2 / ((double)(n - 1));
      S = sqrt(S2);	  






















      


      alpha2 = alpha * alpha;
      mu2 = mu * mu;
      for (j = 0; j < n; j++) {
	for (k = 0; k < n; k++) {
	  //      for (j = 1; j < n; j++) {
	  //	for (k = 0; k < j; k++) {
	  xx = (x[j] - x[k]) / S;
	  xx2 = xx * xx;
	  mux = mu * xx;
	  cosmux = cos(mux);
	  sinmux = sin(mux);
	  z = Faddeeva::erf((C(0.0, 1.0) * xx + alpha2 * mu) / (sqrt(2.0) * alpha), 0.0);
	  term1 = term1 + exp(- xx2 / (2.0 * alpha2)) * (cosmux * (1.0 + creal(z)) - sinmux * cimag(z));
	}
      }
      term1 = term1 / (2.0 * R_pow((double)n, 2.0));
	//      term1 = 1.0 / (double)n + term1 / R_pow((double)n, 2.0);
	  
      for (j = 0; j < n; j++) {
	y = (x[j] - xbar) / S;
	y2 = y * y;
	tmp1 = mu * y / (1.0 + 1.0 / alpha2);
	tmp2 = exp(- y2 / (2.0 * (alpha2 + 1.0)));
	sinmuy = sin(tmp1);
	cosmuy = cos(tmp1);
	z = Faddeeva::erf((alpha2 * mu + C(0.0, 1.0) * y) / (sqrt(2.0) * sqrt(alpha2 + 1.0)), 0.0);

	term2 = term2 + tmp2 * (cosmuy * (1.0 + creal(z)) - sinmuy * cimag(z));
      }
      term2 = R_pow(1.0 + 1.0 / alpha2, -0.5) * exp(-mu2 / (2.0 + 2.0 / alpha2)) * term2 / (double)n;
	  
      term3 = R_pow(1.0 + 2.0 / alpha2, -0.5) * exp(-mu2 / (1.0 + 2.0 / alpha2)) * Rf_pnorm5(alpha * mu / sqrt(1.0 + 2.0 / alpha2), 0.0, 1.0, 1, 0); 

      calphamu = Rf_pnorm5(alpha * mu, 0.0, 1.0, 1, 0);

      statTEP = ((double)n) * (term1 - term2 + term3) / calphamu; 
      



      statistic[0] = statTEP; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
        #include "pvalues/pvalue90.cpp"
      }

// We take the decision to reject or not to reject the null hypothesis H0
      for (i = 0; i < nblevel[0];i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0;   // less
	} else {
	  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
	}
      }

 // If applicable, we free the unused array of pointers
    }

// We return
    return;
          
  }
  
}
