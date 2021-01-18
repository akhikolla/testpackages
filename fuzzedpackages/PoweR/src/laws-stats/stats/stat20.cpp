// Title: Statistique de test de Cabana-Cabana T_{K,l}
// Ref. (book or article): Cabana, A. and Cabana, E. (1994), Goodness-of-Fit and Comparison Tests of the Kolmogorov-Smirnov Type for Bivariate Populations, 
//						   The Annals of Statistics, Vol. 22, No. 3, pp. 1447-1459.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat20(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{K,5}$";
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
//    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double dnorm4(double x, double mu, double sigma, int give_log);
    double *z, *H0, *H1, *H2, *H3, *H4, *H5, *H6, *H7, *H8, *vectoraux2;
    z = new double[n];
    H0 = new double[n];
    H1 = new double[n];
    H2 = new double[n];
    H3 = new double[n];
    H4 = new double[n];
    H5 = new double[n];
    H6 = new double[n];
    H7 = new double[n];
    H8 = new double[n];
    vectoraux2 = new double[n];
    
	double varX=0.0, meanX=0.0, sdX, statTKl, tmp, H3tilde=0.0, H4tilde=0.0, H5tilde=0.0, H6tilde=0.0, H7tilde=0.0, H8tilde=0.0;
    

    for (i=0;i<=(n-1);i++) meanX = meanX + x[i];
    meanX =meanX/(double)n;
    for (i=0;i<=(n-1);i++) varX = varX + R_pow(x[i],2.0);
    varX = ((double)n)*(varX/(double)n - R_pow(meanX,2.0))/(double)(n-1); 
    sdX = sqrt(varX);
    for (i=0;i<=(n-1);i++) z[i] = (x[i]-meanX)/sdX;


    for (i=0;i<=(n-1);i++) {

      H0[i] = 1;
      H1[i] = z[i];
      H2[i] = (R_pow(z[i],2.0) - 1.0)/sqrt(2.0);
      H3[i] = (R_pow(z[i],3.0) - 3.0*z[i])/sqrt(6.0);
      H4[i] = (R_pow(z[i],4.0) - 6.0*R_pow(z[i],2.0) + 3.0)/(2.0*sqrt(6.0));
      H5[i] = (R_pow(z[i],5.0) - 10.0*R_pow(z[i],3.0) + 15.0*z[i])/(2.0*sqrt(30.0));
      H6[i] = (R_pow(z[i],6.0) - 15.0*R_pow(z[i],4.0) + 45.0*R_pow(z[i],2.0) - 15.0)/(12.0*sqrt(5.0));
      H7[i] = (R_pow(z[i],7.0) - 21.0*R_pow(z[i],5.0) + 105.0*R_pow(z[i],3.0) - 105.0*z[i])/(12.0*sqrt(35.0));
      H8[i] = (R_pow(z[i],8.0) - 28.0*R_pow(z[i],6.0) + 210.0*R_pow(z[i],4.0) - 420.0*R_pow(z[i],2.0) + 105.0)/(24.0*sqrt(70.0));

      H3tilde = H3tilde + H3[i];   
      H4tilde = H4tilde + H4[i];  
      H5tilde = H5tilde + H5[i];  
      H6tilde = H6tilde + H6[i];  
      H7tilde = H7tilde + H7[i];  
      H8tilde = H8tilde + H8[i]; 

    }
    
    H3tilde = H3tilde/sqrt((double)n);
    H4tilde = H4tilde/sqrt((double)n);
    H5tilde = H5tilde/sqrt((double)n);
    H6tilde = H6tilde/sqrt((double)n);
    H7tilde = H7tilde/sqrt((double)n);
    H8tilde = H8tilde/sqrt((double)n);
    
    for (i=0;i<=(n-1);i++) {

      vectoraux2[i] = (sqrt(2.0/1.0)*H0[i] + H2[i])*H5tilde + (sqrt(3.0/2.0)*H1[i] + H3[i])*H6tilde + (sqrt(4.0/3.0)*H2[i] + H4[i])*H7tilde + (sqrt(5.0/4.0)*H3[i] + H5[i])*H8tilde + (sqrt(5.0/4.0)*H3[i] + H5[i])*H8tilde;

    }
    
    statTKl = fabs(-dnorm4(z[0],0.0,1.0,0)*H3tilde +(Rf_pnorm5(z[0],0.0,1.0,1,0)-z[0]*dnorm4(z[0],0.0,1.0,0))*H4tilde - dnorm4(z[0],0.0,1.0,0)*vectoraux2[0]);
    for (i=1;i<=(n-1);i++) {
      tmp = fabs(-dnorm4(z[i],0.0,1.0,0)*H3tilde +(Rf_pnorm5(z[i],0.0,1.0,1,0)-z[i]*dnorm4(z[i],0.0,1.0,0))*H4tilde - dnorm4(z[i],0.0,1.0,0)*vectoraux2[i]);
      if (statTKl < tmp) statTKl = tmp;
    }
    statistic[0] = statTKl; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue20.cpp"
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
    delete[] z;
    delete[] vectoraux2;
    delete[] H0;
    delete[] H1;
    delete[] H2;
    delete[] H3;
    delete[] H4;
    delete[] H5;
    delete[] H6;
    delete[] H7;
    delete[] H8;

}

// We return
    return;
   
        
  }
  
}
