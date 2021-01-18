// Title: GEP(t1,t2,t3)
// Ref. (book or article): Alain Desgagné , Pierre Lafaye de Micheaux & Alexandre Leblanc (2013): Test of Normality Against
//             Generalized Exponential Power Alternatives, Communications in Statistics - Theory and Methods, 42:1, 164-190

#include <R.h>
#include "Rmath.h"
typedef void integr_fn(double *x, int n, void *ex);

extern "C" {

  void law34(int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$GEP(t1,t2,t3,crit)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 4;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.5;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 1.e-6;
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
    
// Initialization of the parameters
    double t1, t2, t3, crit;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      t1 = 0.5;
      t2 = 0.0;
      t3 = 1.0;
      crit = 1e-6;
      params[0] = 0.5;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 1e-6;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      t1 = params[0];
      t2 = 0.0;
      t3 = 1.0;
      crit = 1e-6;
      params[1] = 0.0;
      params[2] = 1.0;
      params[3] = 1e-6;
   } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      t1 = params[0];
      t2 = params[1];
      t3 = 1.0;
      crit = 1e-6;
      params[2] = 1.0;
      params[3] = 1e-6;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      t1 = params[0];
      t2 = params[1];
      t3 = params[2];
      crit = 1e-6;
      params[3] = 1e-6;
    } else if  (nbparams[0] == 4) {
      t1 = params[0];
      t2 = params[1];
      t3 = params[2];
      crit = params[3];
    } else {
      error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
//    if (condition) error("?? should not be ?? in lawxx!\n");

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();   
    void dgep(double *y, int n, void *ex);
    double ctenormgep(double t1, double t2, double t3);
    double calsiggep(double cte, double t1, double t2, double t3);
    double calculez0(double cte, double sigma, double t1, double t2, double t3, double crit);
    double calcsup(double cte, double sigma, double z0, double t1, double t2, double t3);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    double Rf_runif(double a, double b);
    double Rf_rbinom(double nin, double pp);
    double *ex, *a, *b, *epsabs, *epsrel, *result, *abserr, *work;
    int *neval, *ier, *limit, *lenw, *last, *iwork;
    ex= new double[5];
    ex[0] = 0.0;
    ex[1] = 0.0;
    ex[2] = t1;
    ex[3] = t2;
    ex[4] = t3;
    a = new double[1];
    b = new double[1];
    epsabs = new double[1];
    epsrel = new double[1];
    result = new double[1];
    abserr = new double[1];
    neval = new int[1];
    ier = new int[1];
    limit = new int[1];
    lenw = new int[1];
    last = new int[1];
   
    epsabs[0] = 0.0001220703;
    epsrel[0] = 0.0001220703;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork = new int[limit[0]];
    work = new double[lenw[0]];
    
    double *wavant, *ww, *u, *z, *zavant, *bb, *w;
    wavant = new double[1];
    ww = new double[n];
    u = new double[n];
    z = new double[n];
    zavant = new double[n];
    bb = new double[n];
    w = new double[1];

    double cte, sigma, gsupp, pourc, z0, extr;
    int ncum=0, lw=0, nsim, jj;
    

   
    cte = ctenormgep(t1, t2, t3);
    
    sigma = calsiggep(cte, t1, t2, t3);
    

    ex[0] = cte;
    ex[1] = sigma;

    z0 =  calculez0(cte, sigma, t1, t2, t3, crit);
    
    
    gsupp = calcsup(cte, sigma, z0, t1, t2, t3);
    
   
    for (i = 0; i < n; i++) z[i] = Rf_runif(0.0, 1.0) * 2.0 * z0 - z0;
    
    for (i = 0; i < n; i++) u[i] = Rf_runif(0.0, 1.0);
    
    
    
    for (i = 0; i < n; i++) zavant[i] = z[i];
    
    dgep(z, n, ex);
    for (i = 0; i < n; i++) bb[i] = z[i] / gsupp;
    
    for (i = 0; i < n ; i++) z[i] = zavant[i];
    
    for (i = 0; i < n; i++) {if (u[i] < bb[i]) lw = lw + 1;}
    
    if (lw > 0) {
      delete[] w;
      w = new double[lw];
      jj = 0;
      for (i = 0; i < n; i++) {if (u[i] < bb[i]) {w[jj] = z[i]; jj = jj + 1;}}
      for (i = 0; i < lw; i++) ww[ncum+i] = w[i];
    }
    
    ncum = ncum + lw;
    
    
    if (0.1 > (double)ncum / (double)n) pourc = 0.1; else pourc = (double)ncum / (double)n;
    
    
    while (ncum < n) {
      nsim = (int)floor((double)(n - ncum) / pourc); 
      
      delete[] z;
      z = new double[nsim];
      for (i = 0; i < nsim; i++) z[i] = Rf_runif(0.0, 1.0) * 2.0 * z0 - z0;
      
      delete[] u;
      u = new double[nsim];
      for (i = 0; i < nsim; i++) u[i] = Rf_runif(0.0, 1.0);
      
      delete[] zavant;
      zavant = new double[nsim];
      for (i = 0; i < nsim; i++) zavant[i] = z[i];
      
      dgep(z, nsim, ex);
      delete[] bb;
      bb = new double[nsim];
      for (i = 0; i < nsim; i++) bb[i] = z[i] / gsupp;
      
      for (i = 0; i < nsim; i++) z[i] = zavant[i];
      
      lw = 0;
      for (i = 0; i < nsim; i++) {if (u[i] < bb[i]) lw = lw + 1;}
      
      if (lw > 0) {
	delete[] w; 
	w = new double[lw];
	jj = 0;
	for (i = 0; i < nsim;i++) {if (u[i] < bb[i]) {w[jj] = z[i]; jj = jj + 1;}}
      }
      
      if (lw > n-ncum) {
	delete[] wavant;
	wavant = new double[lw];
	for (i = 0;i < lw; i++) wavant[i] = w[i];
	lw = n - ncum;
	delete[] w;
	w = new double[lw];
	for (i = 0; i < lw; i++) w[i] = wavant[i];
      }
      
      if (lw > 0) {for (i = 0; i < lw; i++) ww[ncum + i] = w[i];}
      
      ncum = ncum + lw;
    }
    
    
    a[0] = -z0;
    b[0] = z0;
    
    Rdqags(dgep, ex, a, b, epsabs, epsrel, result, abserr, neval, ier,limit, lenw, last,iwork, work);
    extr = 1.0 - result[0];
    
    for (i = 0; i < n; i++) {
      if (Rf_runif(0.0, 1.0) < extr) ww[i] = z0 * (Rf_rbinom(1.0, 0.5) * 2.0 - 1.0); 
    }
    
    
    for (i = 0; i < n; i++) x[i] = ww[i];

    if (setseed[0] == 1) PutRNGstate();
    
// If applicable, we free the unused array of pointers. Then we return.
    delete[] ex;
    delete[] a;
    delete[] b;
    delete[] epsabs;
    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;
    delete[] wavant;
    delete[] ww;
    delete[] u;
    delete[] z;
    delete[] zavant;
    delete[] bb;
    delete[] w;

// We return
    return;
    
  }
 


  void dgep(double *y, int n, void *ex) {
    double cte, sigma, t1, t2, t3, *z;
    int i;
    z = new double[n];
    cte = ((double*)ex)[0];
    sigma = ((double*)ex)[1];
    t1 = ((double*)ex)[2];
    t2 = ((double*)ex)[3];
    t3 = ((double*)ex)[4];
    for (i=0;i<n;i++) z[i] = fabs(y[i]/sigma);
    for (i=0;i<n;i++) {
      y[i] = cte*exp(-0.5*R_pow(z[i],t1))*R_pow(1.0+z[i],-t2)*R_pow(log(exp(1.0)+z[i]),-t3)/sigma;
    }

    delete[] z;
  }

  double ctenormgep(double t1, double t2, double t3) {

    void dgep(double *y, int n, void *ex);
    void Rdqagi(integr_fn f, void *ex, double *bound, int *inf,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    double *ex, *bound, *epsabs, *epsrel, *result, *abserr, *work;
    int *inf, *neval, *ier, *limit, *lenw, *last, *iwork;
    double res;
    ex = new double[5];
    ex[0] = 1.0;
    ex[1] = 1.0;
    ex[2] = t1;
    ex[3] = t2;
    ex[4] = t3;
    bound = new double[1];
    bound[0] = 0.0; //has no meaning since interval is doubly-infinite
    inf = new int[1];
    inf[0] = 2; //interval is doubly-infinite
    epsabs = new double[1];
    epsabs[0] = 0.0001220703;
    epsrel = new double[1];
    epsrel[0] = 0.0001220703;
    result = new double[1];
    abserr = new double[1];
    neval = new int[1];
    ier = new int[1];
    limit = new int[1];
    limit[0] = 100;
    lenw = new int[1];
    lenw[0] = 4 * limit[0];
    last = new int[1];
    iwork = new int[limit[0]];
    work = new double[lenw[0]];

    Rdqagi(dgep, ex, bound, inf, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);

    res = 1.0/result[0];

    delete[] ex;
    delete[] bound;
    delete[] inf;
    delete[] epsabs;
    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;

    return(res);
      
}




  void g(double *y, int n, void *ex) {
    
    void dgep(double *y, int n, void *ex);
    int i;
    double *yavant;
    yavant = new double[n];
    for (i=0;i<n;i++) yavant[i] = y[i];

    dgep(yavant,n,ex);

    for (i=0;i<n;i++) y[i] = yavant[i];

    delete[] yavant;

  }

  double calsiggep(double cte, double t1, double t2, double t3) {

    void g(double *y, int n, void *ex);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);

    double crit = 0.30;

    double *ex, *a, *b, *epsabs, *epsrel, *result, *abserr, *work;
    int *neval, *ier, *limit, *lenw, *last, *iwork;
    double res;
    double sigma = 0.1;
    int flag = 1;
    ex = new double[5];
    ex[0] = cte;
    ex[1] = 1.0;
    ex[2] = t1;
    ex[3] = t2;
    ex[4] = t3;
    a = new double[1];
    b = new double[1];
    a[0] = 0.0;
    epsabs = new double[1];
    epsabs[0] = 0.0001220703;
    epsrel = new double[1];
    epsrel[0] = 0.0001220703;
    result = new double[1];
    abserr = new double[1];
    neval = new int[1];
    ier = new int[1];
    limit = new int[1];
    limit[0] = 100;
    lenw = new int[1];
    lenw[0] = 4 * limit[0];
    last = new int[1];
    iwork = new int[limit[0]];
    work = new double[lenw[0]];


    while (flag==1)
      {
	b[0] = sigma;
	Rdqags(g, ex, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	if (result[0]>crit) flag = 0; else sigma = sigma*1.1;
      }

    res = 1.0/sigma;

    delete[] ex;
    delete[] a;
    delete[] b;
    delete[] epsabs;
    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;

    return(res);

  }



  double inte(double aa, double cte, double sigma, double t1, double t2, double t3) {

    void dgep(double *y, int n, void *ex);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);


    double *ex, *a, *b, *epsabs, *epsrel, *result, *abserr, *work;
    int *neval, *ier, *limit, *lenw, *last, *iwork;
    double res;
    ex = new double[5];
    ex[0] = cte;
    ex[1] = sigma;
    ex[2] = t1;
    ex[3] = t2;
    ex[4] = t3;
    a = new double[1];
    a[0] = 0.0; 
    b = new double[1];
    b[0] = aa; 
    epsabs = new double[1];
    epsabs[0] = 0.0001220703;
    epsrel = new double[1];
    epsrel[0] = 0.0001220703;
    result = new double[1];
    abserr = new double[1];
    neval = new int[1];
    ier = new int[1];
    limit = new int[1];
    limit[0] = 100;
    lenw = new int[1];
    lenw[0] = 4 * limit[0];
    last = new int[1];
    iwork = new int[limit[0]];
    work = new double[lenw[0]];
    
    Rdqags(dgep, ex, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);

    res = result[0];

    delete[] ex;
    delete[] a;
    delete[] b;
    delete[] epsabs;
    delete[] epsrel;
    delete[] result;
    delete[] abserr;
    delete[] neval;
    delete[] ier;
    delete[] limit;
    delete[] lenw;
    delete[] last;
    delete[] iwork;
    delete[] work;

    return(res);

  }

  double calculez0(double cte, double sigma, double t1, double t2, double t3, double crit) {

    double inte(double aa, double cte, double sigma, double t1, double t2, double t3);

    double z0, f;
    int flag;

    z0 = 1.0;
    flag = 1;

    while (flag == 1) {
      f = inte(z0,cte, sigma, t1, t2, t3)*2.0;
      if (f > (1-crit)) flag = 0; else z0 = z0*1.1;
    }
    
    return(z0);
  }



  double calcsup(double cte, double sigma, double z0, double t1, double t2, double t3) {

    void dgep(double *y, int n, void *ex);

    double pas, *range, *ff, *ex, sup, *resultat, res;
    int i, flag, pos, lnrange;

    pas = z0/100.0;
    ex = new double[5];
    ex[0] = cte;
    ex[1] = sigma;
    ex[2] = t1;
    ex[3] = t2;
    ex[4] = t3;
      
    lnrange = 101;
    range = new double[lnrange];
    ff = new double[lnrange];
    resultat = new double[1];

    range[0] = 0.0;
    for (i=1;i<lnrange;i++) range[i] = range[i-1] + pas; 

    flag = 1;
    while (flag == 1) {

      for (i=0;i<lnrange;i++) ff[i] = range[i];
      dgep(ff,lnrange,ex);

      sup = ff[0];
      pos = 0;
      for (i=1;i<lnrange;i++) {

	if (ff[i]>sup) {
	  sup = ff[i];
	  pos = i;
	}

      }

      sup = range[pos];

      if (pas < 0.000001) flag = 0;

      if ((sup-2.0*pas) >= 0.0) {
	lnrange = 401;
	delete[] range;
	range = new double[lnrange];
	range[0] = sup-2.0*pas;
	for (i=1;i<lnrange;i++) range[i] = range[i-1] + pas/100.0; 
      } else {
	lnrange = int((sup+2.0*pas)*100.0/pas)+1;
	delete[] range;
	range = new double[lnrange];
	range[0] = 0.0;
	for (i=1;i<lnrange;i++) range[i] = range[i-1] + pas/100.0; 
      }
      
      delete[] ff;
      ff = new double[lnrange];
      pas = pas/100.0;
    }

    resultat[0] = sup;
    dgep(resultat,1,ex);

    res = resultat[0];

    delete[] range;
    delete[] ff;
    delete[] ex;
    delete[] resultat;

    return res;
 

  }

 
}



