// Title: Statistique de test de Cramer-von Mises pour une distribution de Laplace
// Ref. (book or article): Yen, Vincent C. and Moore, Albert H. (1988) 'Modified goodness-of-fit test for the laplace distribution', 
// Communications in Statistics - Simulation and Computation, 17:1, 275-281.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat43(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$W^2$";
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
      void R_rsort (double* x, int n);
      double plaplace(double y);
      
      double *CvM;
      CvM = new double[n];
      double statCvM, tmp = 0.0, bhat, muhat, sumCvM = 0.0;

      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^
      R_rsort(x, n); 		// we sort the data
      if (n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[n / 2 - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[n / 2];
      }
      
      // calculate b^
      for (i = 0; i < n; i++) {
	tmp = tmp + fabs(x[i] - muhat);
      }
      bhat = tmp / (double)n;
      
      // generate vector CvM
      for (i = 0; i < n; i++) {
	CvM[i] = plaplace((x[i] - muhat) / bhat);
      }
	
      // calculate statCvM
      for (i = 1; i <= n; i++) {
	sumCvM = sumCvM + R_pow((double)(2 * i - 1) / (double)(2 * n) - CvM[i - 1], 2.0);
      }
      statCvM = 1.0 / (double)(12 * n) + sumCvM;

      statistic[0] = statCvM; // Here is the test statistic value

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue43.cpp"
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
      delete[] CvM;
	
    }

    // We return
    return;
        
  }
  
  // In stat42.cpp, we already defined this function so no need to include it here.
  // The cumulative Laplace distribution function with \mu = 0 and \theta = 1
  // fabs = returns the absolute value of x (a negative value becomes positive, positive value is unchanged). 
  // double plaplace(double x) {
    // double temp = 0.5 * exp(-fabs(x));
    // return (x <= 0.0) ? temp : 1.0 - temp;
  // }
  
  
  
// These functions theta(), rho(), imhoffunc(), f(), probQsupx() are copied from package CompQuadForm 
// They are used to computed p-value for Cramer-von Mises, Watson and Anderson-Darling test for Laplace distribution  

  
// Description:
// Distribution function (survival function in fact) of quadratic forms in normal variables using Imhof's method.

  double theta(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2)
    
  { // See Imhof (1961), p.423
    int i, m;
    double sum = 0.0;
    m = lambdalen[0];
    for (i=1;i<=m;i=i+1) sum = sum + h[i-1]*atan(lambda[i-1]*u[0]) + delta2[i-1]*lambda[i-1]*u[0]/(1.0+R_pow(lambda[i-1]*u[0],2.0));
    sum = 0.5 * sum - 0.5*x[0]*u[0];
    return(sum);
  }

  double rho(double *u, double *lambda, int *lambdalen, double *h, double *delta2)
    
  { // See Imhof (1961), p.423
    int i, m;
    double prod = 1.0;
    m = lambdalen[0];
    for (i=1;i<=m;i=i+1)  prod = prod * R_pow(1.0+R_pow(lambda[i-1]*u[0],2.0),0.25*h[i-1])*exp(0.5*delta2[i-1]*R_pow(lambda[i-1]*u[0],2.0)/(1.0 + R_pow(lambda[i-1]*u[0],2.0)));
    return(prod);
  }
  

  double imhoffunc(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2)
    
  { // This is the function under the integral sign in equation (3.2), Imhof (1961), p.422
    double theta(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2);
    double rho(double *u, double *lambda, int *lambdalen, double *h, double *delta2);
    double res;
    res = (sin(theta(u,lambda,lambdalen,h,x,delta2)))/(u[0]*rho(u,lambda,lambdalen,h,delta2));
    return(res);
  }
  
  // Already defined in stat29.cpp
  //  typedef void integr_fn(double *x, int n, void *ex);

  void f(double *x, int n, void *ex) 

  {
    double imhoffunc(double *u, double *lambda, int *lambdalen, double *h, double *x, double *delta2);
    int i;
    double *xx;
    xx = new double[1];
    xx[0] = ((double*)ex)[0];
    int *lambdalen;
    lambdalen = new int[1];
    lambdalen[0] = (int)(((double*)ex)[1]);
    double *lambda;
    lambda = new double[lambdalen[0]];
    for (i=1;i<=lambdalen[0];i=i+1) lambda[i-1] = ((double*)ex)[i+1];
    double *h;
    h = new double[lambdalen[0]];
    for (i=1;i<=lambdalen[0];i=i+1) h[i-1] = ((double*)ex)[lambdalen[0]+i+1];
    double *delta2;
    delta2 = new double[lambdalen[0]];
    for (i=1;i<=lambdalen[0];i=i+1) delta2[i-1] = ((double*)ex)[2*lambdalen[0]+i+1];
    double *u;
    u = new double[1];
    for (i=1;i<=n;i=i+1) {
      u[0] = x[i-1];
      x[i-1] =  imhoffunc(u, lambda, lambdalen, h, xx, delta2);
    }

    delete[] xx;
    delete[] lambdalen;
    delete[] lambda;
    delete[] h;
    delete[] delta2;
    delete[] u;

	// debug
	return;
	// end of debug

  }
  

  void probQsupx(double *x, double *lambda, int *lambdalen, double *h, double *delta2, double *Qx, double *epsabs, double *epsrel, int *limit)
    
  {
    int i;
    void f(double *x, int n, void *ex); 
    void Rdqagi(integr_fn f, void *ex, double *bound, int *inf,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    double *ex;
    ex = new double[2+3*lambdalen[0]]; 
    ex[0] = x[0];
    ex[1] = (double)lambdalen[0];
    for (i=1;i<=lambdalen[0];i=i+1) ex[i+1] = lambda[i-1];
    for (i=1;i<=lambdalen[0];i=i+1) ex[lambdalen[0]+i+1] = h[i-1];
    for (i=1;i<=lambdalen[0];i=i+1) ex[2*lambdalen[0]+i+1] = delta2[i-1];
    double *bound;  
    bound = new double[1];
    bound[0] = 0.0;
    int *inf;  
    inf = new int[1];
    inf[0] = 1;
    // double *epsabs;  
    // epsabs = new double[1];
    // *(epsabs+0) = 0.000001; //0.0001220703;
    // double *epsrel;  
    // epsrel = new double[1];
    // *(epsrel+0) = 0.000001; //0.0001220703;
    double *result;  
    result = new double[1];
    double *abserr;  
    abserr = new double[1];
    int *neval;  
    neval = new int[1];
    int *ier;  
    ier = new int[1];
    // int *limit;  
    //   limit = new int[1];
    // *(limit+0) = 10000;
    int *lenw;  
    lenw = new int[1];
    *(lenw+0) = 4 * *(limit+0);
    int *last;  
    last = new int[1];
    int *iwork;  
    iwork = new int[*(limit+0)];
    double *work;  
    work = new double[*(lenw+0)];
    Rdqagi(f,ex,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,iwork,work);
    Qx[0] = 0.5 + result[0]/M_PI;
    
    epsabs[0] = abserr[0];

    delete[] ex;
    delete[] bound;
    delete[] inf;
    //    delete[] epsabs;
    //    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    //    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;
    
	return;
  }
  
// End of theta(), rho(), imhoffunc(), f(), probQsupx() functions !!!
  
  
}
