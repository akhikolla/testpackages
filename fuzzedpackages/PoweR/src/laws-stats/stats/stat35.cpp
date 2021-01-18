// Title: The Rn test for normality
// Ref. (book or article): Desgagne, A., Lafaye de Micheaux, P. and Leblanc, A. (2013), Test of Normality Against Generalized Exponential Power Alternatives, 
//                         Communications in Statistics - Theory and Methods, 42, 164--190.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat35(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$R_n$";
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
      for (i = j; i < 50; i++) name[i][0] = space[0];
      return;
    }

    if (n > 3) {
// Computation of the value of the test statistic
    double *y;
    y = new double[n];
    double varpopX = 0.0, meanX = 0.0, sdX, r1 = 0.0, r2 = 0.0, r3 = 0.0, Rn;

    for (i = 0; i <= (n - 1); i++) meanX = meanX + x[i];
    meanX = meanX / (double)n;
    for (i = 0; i <= (n - 1); i++) varpopX = varpopX + R_pow(x[i], 2.0);
    varpopX = varpopX / (double)n - R_pow(meanX, 2.0); 
    sdX = sqrt(varpopX);
    for (i = 0; i <= (n - 1); i++) y[i] = (x[i] - meanX) / sdX;

    // Formulas given in our paper p. 169
    for (i = 0; i <= (n - 1); i++) {
      r1 = r1 + R_pow(y[i], 2.0) * log(fabs(y[i]));
      r2 = r2 + log(1.0 + fabs(y[i]));
      r3 = r3 + log(log(2.71828182846 + fabs(y[i])));
    }
    r1 = 0.18240929 - 0.5 * r1 / (double)n;
    r2 = 0.5348223 - r2 / (double)n;
    r3 = 0.20981558 - r3 / (double)n;
 
    // Formula given in our paper p. 170
    Rn = (double)n * ((r1 * 1259.04213344 - r2 * 32040.69569026 + r3 * 85065.77739473) * r1 + (-r1 * 32040.6956903 + r2 * 918649.9005906 - r3 * 2425883.3443201) * r2 + (r1 * 85065.7773947 - r2 * 2425883.3443201 + r3 * 6407749.8211208) * r3);
    
    statistic[0] = Rn; // Here is the test statistic value

    if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue35.cpp"
    }

// We take the decision to reject or not to reject the null hypothesis H0
     for (i = 0; i <= (nblevel[0] - 1); i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
	    if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers
    delete[] y;

}

// We return
    return;
   
        
  }
  
}
