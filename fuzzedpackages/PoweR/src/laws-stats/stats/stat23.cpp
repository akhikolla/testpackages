// Title: Statistique de test de Shapiro-Wilk modifiée par Rahman-Govindarajulu 
// Ref. (book or article): (Voir le package nortest) ?
//							Rahman, M.M. and Govindarajulu, Z. (1997), A modification of the test of Shapiro and Wilk for normality, 
//							Journal of Applied Statistics, Vol. 24, Issue 2, pp. 219-236.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat23(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 4;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$\\tilde{W}$";
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
    //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double qnorm(double p, double mean, double sd, int lower_tail, int log_p);
    double dnorm4(double x, double mu, double sigma, int give_log);
    double *mi, *fi, *xs, *aux1, *aux2, *aux3, *aux4, *aistar, *ai;
    xs = new double[n];
    mi = new double[n];
    fi = new double[n];
    aux1 = new double[n];
    aux2 = new double[n];
    aux3 = new double[n];
    aux4 = new double[n];
    aistar = new double[n];
    ai = new double[n];
    double meanX=0.0, statWRG=0.0, norm2=0.0, aux6=0.0;
 

    for (i=1;i<=n;i++) { 
      mi[i-1] = qnorm((double)i/(double)(n+1),0.0,1.0,1,0);
      fi[i-1] = dnorm4(mi[i-1],0.0,1.0,0);
      aux2[i-1] = 2*mi[i-1]*fi[i-1];
    }

    aux1[0] = 0.0;
    for (i=1;i<=(n-1);i++) { 
      aux1[i] = mi[i-1]*fi[i-1];
    }

    for (i=0;i<=(n-2);i++) { 
      aux3[i] = mi[i+1]*fi[i+1];
    }
    aux3[n-1] = 0.0;

    for (i=0;i<=(n-1);i++) {
      aux4[i] = aux1[i] - aux2[i] +aux3[i];
      aistar[i] = -((double)(n+1)*(n+2))*fi[i]*aux4[i];
      norm2 = norm2 + R_pow(aistar[i],2.0);
    }

    for (i=0;i<=(n-1);i++) ai[i] = aistar[i]/sqrt(norm2);

    for (i=0;i<=(n-1);i++) xs[i] = x[i];
    R_rsort (xs,n); // We sort the data

    for (i=0;i<=(n-1);i++) meanX = meanX + xs[i];
    meanX =meanX/(double)n;

    for (i=0;i<=(n-1);i++) {
      aux6 = aux6 + R_pow(xs[i]-meanX,2.0);
    }
    for (i=0;i<=(n-1);i++) {
      statWRG = statWRG + ai[i]*xs[i];
    }
    statWRG = R_pow(statWRG,2.0)/aux6;

    statistic[0] = statWRG; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue23.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
   for (i=0;i<=(nblevel[0]-1);i++) {   
      if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values
      } else {
	//	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
   
// If applicable, we free the unused array of pointers
   delete[] xs;
   delete[] mi;
   delete[] fi;
   delete[] aux1;
   delete[] aux2;
   delete[] aux3;
   delete[] aux4;
   delete[] aistar;
   delete[] ai;

}

// We return
   return;
   
        
  }
  
}
