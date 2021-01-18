// Title: The X_{APD} test
// Ref. (book or article): A Powerful and Interpretable Alternative to the Jarque-Bera Test of Normality Based on 2nd-Power Skewness and Kurtosis, using the Rao’s score test on the APD family, Journal of Applied Statistics

#include <R.h>
#include "Rmath.h"

extern "C" {
  void stat36(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;    

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$X_{APD}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 0;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).

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

    if (n>3) {
// Computation of the value of the test statistic
    double qchisq(double p, double df, int lower_tail, int log_p);
    double pchisq(double q, double df, int lower_tail, int log_p);
    double *zchap, C = 0.0, S = 0.0, tmp = 0.0;
    zchap = new double[n];
    double gamma = 0.577215664901533, muchap = 0.0, sigchap, sigchap2 = 0.0, statvalue;

    for (i = 0; i <= (n - 1); i++) muchap = muchap + x[i];
    muchap = muchap / (double)n;

    for (i = 1; i <= n; i++) sigchap2 = sigchap2 + R_pow(x[i - 1] - muchap, 2.0);
    sigchap2 = sigchap2 / (double)n;
    sigchap = sqrt(sigchap2);
    for (i = 0; i <= (n - 1); i++) zchap[i] = (x[i] - muchap) / sigchap;

    
    for (i = 0; i <= (n - 1); i++) {
      if (zchap[i] != 0) {
	tmp = R_pow(zchap[i], 2.0);
	C = C + tmp * log(tmp);
	S = S + zchap[i] * fabs(zchap[i]);
      }
    }

    C = C / (double)n; // This is 2K_2
    S = S / (double)n; // This is B_2

    double kurtcond = R_pow(C - 2 * R_pow(S, 2.0), 1.0 / 3.0);
    double skewtranssquare = (double)n * R_pow(S, 2.0) / ((3.0 - 8.0 / M_PI) * (1.0 - 1.9 / ((double)n)));
    double espku = R_pow(2.0 - log(2.0) - gamma, 1.0 / 3.0) * (1.0 - 1.026 / ((double)n));
    double varku = ((3.0 * M_PI * M_PI / 2.0 - 14.0) / R_pow(3.0 * R_pow(2.0 - log(2.0) - gamma, 2.0 / 3.0) , 2.0)) * (1.0 - 2.25 / R_pow((double)n, 0.8)) / ((double)n);
    double kurtcondtranssquare = R_pow(kurtcond - espku, 2.0) / varku;
														      
    statvalue = skewtranssquare + kurtcondtranssquare;

    //    statvalue = ((double)n) * (R_pow(C - (2.0 - log(2.0) - gamma), 2.0) / (1.5 * M_PI * M_PI - 14.0) + S * S / (3.0 - 8.0 / M_PI));

    statistic[0] = statvalue; // Here is the test statistic value

    if (pvalcomp[0] == 1) {
      // If possible, computation of the p-value.
      #include "pvalues/pvalue36.cpp"
    }

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i < nblevel[0]; i++) {
      
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers
    delete[] zchap;

}

// We return
    return;
    
        
  }
 
 
}

