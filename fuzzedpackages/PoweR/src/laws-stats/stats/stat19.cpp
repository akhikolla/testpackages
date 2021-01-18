// Title: Statistique de test de Cabana-Cabana T_{S,l}
// Ref. (book or article): Cabana, A. and Cabana, E. (1994), Goodness-of-Fit and Comparison Tests of the Kolmogorov-Smirnov Type for Bivariate Populations, 
//						   The Annals of Statistics, Vol. 22, No. 3, pp. 1447-1459.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat19(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{S,5}$";
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
    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double dnorm4(double x, double mu, double sigma, int give_log);
    double *zdata, *vectoraux1;
    zdata = new double[n];
    vectoraux1 = new double[n];
    double varX=0.0, meanX=0.0, sdX, statTSl, tmp, meanH3=0.0, meanH4=0.0, meanH5=0.0, meanH6=0.0, meanH7=0.0, meanH8=0.0;

    for (i=0;i<=(n-1);i++) meanX = meanX + x[i];
    meanX =meanX/(double)n;
    for (i=0;i<=(n-1);i++) varX = varX + R_pow(x[i],2.0);
    varX = ((double)n)*(varX/(double)n - R_pow(meanX,2.0))/(double)(n-1); 
    sdX = sqrt(varX);
    for (i=0;i<=(n-1);i++) zdata[i] = (x[i]-meanX)/sdX;


    for (i=0;i<=(n-1);i++) {

      meanH3 = meanH3 + (R_pow(zdata[i],3.0) - 3.0*zdata[i])/sqrt(6.0);   

      meanH4 = meanH4 + (R_pow(zdata[i],4.0) - 6.0*R_pow(zdata[i],2.0) + 3.0)/(2.0*sqrt(6.0));  

      meanH5 = meanH5 + (R_pow(zdata[i],5.0) - 10.0*R_pow(zdata[i],3.0) + 15.0*zdata[i])/(2.0*sqrt(30.0));  

      meanH6 = meanH6 + (R_pow(zdata[i],6.0) - 15.0*R_pow(zdata[i],4.0) + 45.0*R_pow(zdata[i],2.0) - 15.0)/(12.0*sqrt(5.0));  

      meanH7 = meanH7 + (R_pow(zdata[i],7.0) - 21.0*R_pow(zdata[i],5.0) + 105.0*R_pow(zdata[i],3.0) - 105.0*zdata[i])/(12.0*sqrt(35.0));  

      meanH8 = meanH8 + (R_pow(zdata[i],8.0) - 28.0*R_pow(zdata[i],6.0) + 210.0*R_pow(zdata[i],4.0) - 420.0*R_pow(zdata[i],2.0) + 105.0)/(24.0*sqrt(70.0)); 

      vectoraux1[i] = 0.0;

    }
    
    meanH3 = meanH3/sqrt((double)n);
    meanH4 = meanH4/sqrt((double)n);
    meanH5 = meanH5/sqrt((double)n);
    meanH6 = meanH6/sqrt((double)n);
    meanH7 = meanH7/sqrt((double)n);
    meanH8 = meanH8/sqrt((double)n);
    
    for (i=0;i<=(n-1);i++) {

      vectoraux1[i] = meanH4 + meanH5*zdata[i]/sqrt(2.0) + meanH6*(R_pow(zdata[i],2.0) - 1.0)/sqrt(6.0) + meanH7*(R_pow(zdata[i],3.0) - 3.0*zdata[i])/(2.0*sqrt(6.0)) + meanH8*(R_pow(zdata[i],4.0) - 6.0*R_pow(zdata[i],2.0) + 3.0)/(2.0*sqrt(30.0));

    }
    
    statTSl = fabs(pnorm(zdata[0],0.0,1.0,1,0)*meanH3 - dnorm4(zdata[0],0.0,1.0,0)*vectoraux1[0]);
    for (i=1;i<=(n-1);i++) {
      tmp = fabs(pnorm(zdata[i],0.0,1.0,1,0)*meanH3 - dnorm4(zdata[i],0.0,1.0,0)*vectoraux1[i]);
      if (statTSl < tmp) statTSl = tmp;
    }
    statistic[0] = statTSl; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue19.cpp"
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
    delete[] zdata;
    delete[] vectoraux1;

}

// We return
    return;
   
        
  }
  
}
