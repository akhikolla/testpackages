// Title: Statistique de test de Coin \beta_3^2
// Ref. (book or article): Coin, D. (2008), A goodness-of-fit test for normality based on polynomial regression, 
//						   Computational Statistics & Data Analysis, Vol. 52, pp. 2185-2198.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat30(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$\\beta_3^2$";
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
      void nscor2(double *s, int *n, int *n2);
      //      double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
      double *z, *sp, *a;
      int *M;
      z = new double[n];
      M = new int[1];
      M[0] = n/2; // integer division
      sp = new double[M[0]];
      a = new double[n];
      double varX=0.0, meanX=0.0, sdX, statbeta32, term1=0.0, term2=0.0, term3=0.0, term4=0.0, term6=0.0;
      
      for (i=0;i<=(n-1);i++) {
	    meanX = meanX + x[i];
	  }
      meanX =meanX/(double)n;
      
	  for (i=0;i<=(n-1);i++) {
	    varX = varX + R_pow(x[i],2.0);
	  }
      varX = ((double)n)*(varX/(double)n - R_pow(meanX,2.0))/(double)(n-1); 
      sdX = sqrt(varX);
      
	  for (i=0;i<=(n-1);i++) {
	    z[i] = (x[i]-meanX)/sdX;
	  }

      R_rsort(z,n);
      
      nscor2(sp,xlen,M);
      
      if ((n%2) == 0) {
	
	// c(-value, rev(value))
	for (i=0;i<=((n/2)-1);i++) a[i] = -sp[i]; 
	for (i=(n/2);i<=(n-1);i++) a[i] = sp[n-i-1];
	
      } else { // n is odd
	
	// c(-value, 0, rev(value))
	for (i=0;i<=((n/2)-1);i++) a[i] = -sp[i]; 
	a[(n/2)] = 0.0;
	for (i=((n/2)+1);i<=(n-1);i++) a[i] = sp[n-i-1];
	
      }
      
      for (i=0;i<=(n-1);i++) {
	
		term1 = term1 + R_pow(a[i],4.0);
		term2 = term2 + a[i]*z[i];
		term3 = term3 + R_pow(a[i],2.0);
		term4 = term4 + R_pow(a[i],3.0)*z[i];
		term6 = term6 + R_pow(a[i],6.0);
	
      }
      
    statbeta32 = R_pow((term1*term2-term3*term4)/(term1*term1-term3*term6),2.0);
      
    statistic[0] = statbeta32; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue30.cpp"
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
      delete[] z;
      delete[] sp;
      delete[] M;
      delete[] a;

    }

// We return
    return;
   
        
  }


static double correc(int, int);

void nscor2(double *s, int *n, int *n2) {

/*   algorithm as 177.3, applied statistics, v.31, 161-165, 1982.

     calculates approximate expected values of normal order statistics.
     claimed accuracy is 0.0001, though usually accurate to 5-6 dec.

 ***  N.B. This routine was NOT in double precision All constants were f ***

     Arguments:

     s(n2)   = output, the first n2 expected values.
     n	     = input, the sample size.
     n2	     = input, the number of order statistics required; must
		      be <= n/2.

	 ier removed REW Mar 2001

     ier   = output, error indicator
		   = 0 if no error detected
		   = 1 if n <= 1.
		   = 2 if n > 2000, in which case the order statistics
			  are still calculated, but may be inaccurate.
		   = 3 if n2 > n/2 (n.b. this differs from the
			  published algorithm which returns an error
			  if n2 is not equal to n/2.)

     Calls qnorm() [from R] which is an improvement of
     ppnd = applied statistics algorithm 111.
     An alternative is ppnd7 in algorithm AS 241.
*/

    /* Initialized data */

    const double
	eps[4] = { .419885,.450536, .456936, .468488 },
	dl1[4] = { .112063,.12177,  .239299, .215159 },
	dl2[4] = { .080122,.111348,-.211867,-.115049 },
	gam[4] = { .474798,.469051, .208597, .259784 },
	lam[4] = { .282765,.304856, .407708, .414093 };
    const double bb = -.283833;
    const double d  = -.106136;
    const double b1 = .5641896;


    /* Local variables */
    int i, k;
    double e1, e2, ai, an;

    /* input parameter checks. */

    if (*n2 > *n / 2) {
		error("\nn2>n");
    }
    if (*n <= 1) {
		error("\nn<=1");
    }
    if (*n > 2000) {
		warning("\nValues may be inaccurate because of the size of N");
    }

    s[0] = b1;
    if (*n == 2) {
	return;
    }

/*	calculate normal tail areas for first 3 order statistics. */

    an = (double) (*n);
    k = 3;
    if (*n2 < k)
	k = *n2;
    /* k := min(3, *n2) */
    for (i = 0; i < k; ++i) {
	ai = (double) i+1;
	e1 = (ai - eps[i]) / (an + gam[i]);
	e2 = pow((double) e1, (double) lam[i]);
	s[i] = e1 + e2 * (dl1[i] + e2 * dl2[i]) / an - correc(i+1, *n);
    }
    if (*n2 > k) {

/*	calculate normal areas for other cases. */

	for (i = 4-1; i < *n2; ++i) {
	    ai = (double) i+1;
	    e1 = (ai - eps[3]) / (an + gam[3]);
	    e2 = pow((double) e1, (double) lam[3] + bb / (ai + d));
	    s[i] = e1 + e2 * (dl1[3] + e2 * dl2[3]) / an - correc(i+1, *n);
	}
    }
/*	convert tail areas to normal deviates. */

    for (i = 0; i < *n2; ++i)
	s[i] = - qnorm(s[i], 0., 1., 1, 0);

    return;
} /* nscor2 */


static double correc(int i, int n)
{
/*	calculates correction for tail area of the i-th largest of n
	order statistics. */

    const double
	c1[7] = { 9.5,28.7,1.9,0.,-7.,-6.2,-1.6 },
	c2[7] = { -6195.,-9569.,-6728.,-17614.,-8278.,-3570., 1075. },
	c3[7] = { 93380.,175160.,410400.,2157600.,2.376e6,2.065e6,2.065e6 };
    const double mic = 1e-6;
    const double c14 = 1.9e-5;

    double an;

    if (i * n == 4)		return c14;
    if (i < 1 || i > 7)		return 0;
    if (i != 4 && n > 20)	return 0;
    if (i == 4 && n > 40)	return 0;
    /* else : */
    an = (double) n;
    an = 1. / (an * an);
    i--;
    return((c1[i] + an * (c2[i] + an * c3[i])) * mic);
} /* correc */



  
}
