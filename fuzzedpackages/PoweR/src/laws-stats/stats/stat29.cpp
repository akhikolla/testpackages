// Title: Statistique de test de del Barrio-Cuesta-Albertos-Matran-Rodriguez-Rodriguez quantile
// Ref. (book or article): Barrio, E. del, Cuesta-Albertos, J., Matran, C. and Rodriguez-Rodriguez, J. (1999), 
//						   Tests of goodness-of-fit based on the L_2-Wasserstein distance,
//						   The Annals of Statistics, Vol. 27, pp. 1230-1239. 

#include <R.h>
#include "Rmath.h"
typedef void integr_fn(double *x, int n, void *ex);

extern "C" {
 
  void stat29(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$BCMR$";
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
      void f29(double *z, int n, void *ex);
      void Rdqags(integr_fn f29, void *ex, double *a, double *b,
	              double *epsabs, double *epsrel,
	              double *result, double *abserr, int *neval, int *ier,
	              int *limit, int *lenw, int *last,
	              int *iwork, double *work);
      void R_rsort (double* x, int n);
      double pchisq(double q, double df, int lower_tail, int log_p);
      double *xtmp;
      xtmp = new double[n];
      double *a, *b, *epsabs, *epsrel, *result, *ex, *abserr, *work;
      int *last, *limit, *lenw, *ier, *neval, *iwork;
      epsabs = new double[1];
      epsrel = new double[1];
      ex = new double[0];
      limit = new int[1];
      lenw = new int[1];
      a = new double[1];
      b = new double[1];
      result = new double[1];
      abserr = new double[1];
      last = new int[1];
      ier = new int[1];
      neval = new int[1];
        
      epsabs[0] = 0.0001220703;
      epsrel[0] = 0.0001220703;
      
      limit[0] = 100;
      lenw[0] = 4 * limit[0];
      
      iwork = new int[limit[0]];
      work = new double[lenw[0]];
      
      double statBCMR=0.0, m2=0.0, meanX=0.0, val;

      for (i=0;i<=(n-1);i++) meanX = meanX + x[i];
      meanX = meanX/(double)n;
      for (i=0;i<=(n-1);i++) {
		xtmp[i] = x[i];
		m2 = m2 + R_pow(x[i]-meanX,2.0);
      }
      m2 = m2/(double)n;

      R_rsort(xtmp,n); // We sort the data

      
      for (i=1;i<=n;i++) {

		a[0] = (double)(i-1)/(double)n;
		b[0] = (double)i/(double)n;

		Rdqags(f29, ex, a, b, epsabs, epsrel, result, abserr, neval, ier,limit, lenw, last,iwork, work);
		val = result[0];

		statBCMR = statBCMR + xtmp[i-1]*val;

      }


      statBCMR = 1.0 - R_pow(statBCMR,2.0)/m2;
      
      // Here is the test statistic value
      statistic[0] = statBCMR; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue29.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
	//	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }

 // If applicable, we free the unused array of pointers
    delete[] xtmp;
    delete[] a;
    delete[] b;
    delete[] epsabs;
    delete[] epsrel;
    delete[] result;
    delete[] ex;
    delete[] abserr;
    delete[] work;
    delete[] last;
    delete[] limit;
    delete[] lenw;
    delete[] ier;
    delete[] neval;
    delete[] iwork;

    }
    
// We return
    return;
   
    
  }

 
  void f29(double *z, int n, void *ex) {
    double qnorm(double q, double mean, double sd, int lower_tail, int log_p);
    int i;    
    for (i=0;i<=(n-1);i++) z[i] = qnorm(z[i],0.0,1.0,1,0);
  }

  
}
