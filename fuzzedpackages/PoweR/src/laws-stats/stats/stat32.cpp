// Title: Statistique de test de Martinez-Iglewicz
// Ref. (book or article): Martinez, J. and Iglewicz, B. (1981), \emph{A test for departure from normality based on a biweight estimator of scale}, 
//						   Biometrika, Vol. 68, No. 1, pp. 331-333.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat32(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$I_n$";
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
      double pchisq(double q, double df, int lower_tail, int log_p);
      double *xtmp, *aux1, *z, z2;
      xtmp = new double[n];
      aux1 = new double[n];
      z = new double[n];
      double statIn, M, A, Sb2, term1=0.0, term2=0.0, term3=0.0;


      for (i=0;i<n;i++) {
	    xtmp[i] = x[i];
	  }
      
      R_rsort(xtmp,n); // We sort the data
      if ((n%2) == 0) {
	    M = (xtmp[n/2]+xtmp[n/2-1])/2.0; 
	  } else {
	    M = xtmp[n/2]; // sample median
	  }

      for (i=0;i<n;i++) {
	aux1[i] = x[i] - M;
	xtmp[i] = fabs(aux1[i]);
      }

      R_rsort(xtmp,n); // We sort the data
      if ((n%2) == 0) {
	A = (xtmp[n/2]+xtmp[n/2-1])/2.0; 
      } else {
	A = xtmp[n/2]; // sample median
      }

      A = 9.0 * A;

      for (i=0;i<n;i++) {
	z[i] = aux1[i] / A;
	if (fabs(z[i]) < 1.0) {
	  z2 = R_pow(z[i], 2.0);
	  term1 = term1 + R_pow(aux1[i], 2.0) * R_pow(1.0 - z2, 4.0);
	  term2 = term2 + (1.0 - z2) * (1.0 - 5.0 * z2);		
	}
	term3 = term3 + R_pow(aux1[i], 2.0);
      }
      // Problem here!! term2 might be equal to 0 (implying a division by 0 afterwards)
      // The test statistic is not well defined in this paper (neither in others I looked at)
      // BUT
      /*
J'ai regardé cet estimateur. Les zi, si k=1, sont des données centrées réduites. Avec ce petit code, ça permet de voir comment le z se comporte. Ici je pose k=1 et la moitié de l'échantillon (pour n pair) a |z|<1. Si on pose k=9 comme ils font, alors c'est presque la totalité de l'échantillon qui a |z|<1. Le cas où tous les |z| sont >= 1 est impossible je crois.

n <- 8
k <- 1
a <- NULL
for (i in 1:100){
  x <- rnorm(n)
  z1 <- x-median(x);
  z <- z1/median(abs(z1));
  a <- c(a,sum(abs(z/k)<1));
}
sort(a)
       */

      Sb2 = ((double)n) * term1 / R_pow(term2, 2.0);

      statIn = (term3 / ((double)(n - 1))) / Sb2;
      
      statistic[0] = statIn; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue32.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
	  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }

 // If applicable, we free the unused array of pointers
      delete[] xtmp;
      delete[] aux1;
      delete[] z;
    }
    
// We return
    return;
   
    
  }
  
}
