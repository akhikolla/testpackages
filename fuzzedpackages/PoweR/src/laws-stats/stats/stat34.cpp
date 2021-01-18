// Title: Statistique de test de Zhang Q*
// Ref. (book or article): Zhang, P (1999), Omnibus test of normality using the Q statistic, 
//						   Journal of Applied Statistics, Vol. 26, Issue 4, pp. 519-528. 

// A VOIR !!! ON N'A PAS LA P-VALEUR!!

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat34(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
 // Here, INDICATE the name of your statistic
      const char *nom = "Q*";
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
    void R_rsort (double* x, int n);
    double qnorm(double p, double mean, double sd, int lower_tail, int log_p);
    //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double *u, *xs, *a, *b;
    u = new double[n];
    a = new double[n];
    b = new double[n];
    xs = new double[n];
    double Qstar, q1star=0.0, q2star=0.0, term=0.0;
 

    for (i=0;i<=(n-1);i++) *(xs+i) = *(x+i);
    R_rsort (xs,n); // We sort the data

    for (i=1;i<=n;i++) u[i-1] = qnorm(((double)i-0.375)/((double)n+0.25),0.0,1.0,1,0);

    for (i=2;i<=n;i++) {
      a[i-1] = 1.0/(((double)(n-1))*(u[i-1]-u[0]));
      term = term + a[i-1];
    }
    a[0] = -term;
    // ATTENTION!! Il faut que n>8 ???
    b[0] = 1.0/(((double)(n-4))*(u[0]-u[4]));
    b[n-1] = -b[0];
    b[1] = 1.0/(((double)(n-4))*(u[1]-u[5]));
    b[n-2] = -b[1];
    b[2] = 1.0/(((double)(n-4))*(u[2]-u[6]));
    b[n-3] = -b[2];
    b[3] = 1.0/(((double)(n-4))*(u[3]-u[7]));
    b[n-4] = -b[3];
    for (i=5;i<=(n-4);i++) {
      b[i-1] = (1.0/(u[i-1]-u[i+3]) - 1.0/(u[i-5]-u[i-1]))/((double)(n-4));
    }

    for (i=1;i<=n;i++) { 
      q1star = q1star - a[i-1]*xs[n-i];
      q2star = q2star - b[i-1]*xs[n-i];
    }

    
    Qstar= log(q1star/q2star);

    *(statistic+0) = Qstar; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue34.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=1;i<=nblevel[0];i++) {
     
    if (usecrit[0] == 1) { // We use the provided critical values
	if (alter[0] == 0) { if (statistic[0] > critvalR[i-1] || statistic[0] < critvalL[i-1]) decision[i-1] = 1; else decision[i-1] = 0; // two-sided
	} else if (alter[0] == 1) {  if (statistic[0] < critvalL[i-1]) decision[i-1] = 1; else decision[i-1] = 0;  // less
	} else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i-1]) decision[i-1] = 1; else decision[i-1] = 0; } } // greater
      } else {
	//     if (pvalue[0] < level[i-1]) decision[i-1] = 1; else decision[i-1] = 0; // We use the p-value
      }
    }

// If applicable, we free the unused array of pointers
  delete[] u;
  delete[] a;
  delete[] b;
  delete[] xs;

}

// We return
  return;
   
        
  }
  
}
