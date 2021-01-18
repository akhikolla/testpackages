// Title: Statistique de test de D'Agostino 
// Ref. (book or article): An omnibus test of normality for moderate and large size samples. Biometrika 58 (1971), 341--348.

// A VOIR !!! ON N'A PAS LA P-VALEUR!!

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat24(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
 // Here, INDICATE the name of your statistic
      const char *nom = "$D$";
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
    double *xs;
    xs = new double[n];
    double varX=0.0, meanX=0.0, T=0.0, D, statDa, valcritleft=0.0, valcritright=0.0;
 
    for (i=0;i<=(n-1);i++) xs[i] = x[i];
    R_rsort (xs,n); // We sort the data
    for (i=0;i<=(n-1);i++) meanX = meanX + xs[i];
    meanX =meanX/(double)n;
    for (i=0;i<=(n-1);i++) varX = varX + R_pow(xs[i],2.0);
    varX = varX/(double)n -R_pow(meanX,2.0); // !! var. pop. ici
    for (i=1;i<=n;i++) T = T + ((double)i-0.5*((double)(n+1)))*xs[i-1];
    D = T/(((double)(n*n))*sqrt(varX));
    statDa = sqrt((double)n)*(D - 0.28209479)/0.02998598;

    statistic[0] = statDa; // Here is the test statistic value


if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue24.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
      for (i=0;i<=(nblevel[0]-1);i++) {
      if (level[i] == 0.05) { // calculés avec 6.10^7 boucles de Monte-Carlo
	
	if (n == 20) {valcritleft = -3.044; valcritright = 0.634;}
	if (n == 50) {valcritleft = -2.740; valcritright = 1.058;}
	if (n == 100) {valcritleft = -2.543; valcritright = 1.312;}
	if (n == 200) {valcritleft = -2.387; valcritright = 1.500;}
	if (n == 500) {valcritleft = -2.238; valcritright = 1.669;}
	
	if ( (statDa <= valcritleft) || (statDa >= valcritright)) decision[i] = 1; else decision[i] = 0; 
	
      }
      
      if (level[i] == 0.1) { // calculés avec 6.10^7 boucles de Monte-Carlo
	
	if (n == 20) {valcritleft = -2.439; valcritright = 0.563;}
	if (n == 50) {valcritleft = -2.212; valcritright = 0.937;}
	if (n == 100) {valcritleft = -2.070; valcritright = 1.144;}
	if (n == 200) {valcritleft = -1.958; valcritright = 1.293;}
	if (n == 500) {valcritleft = -1.850; valcritright = 1.424;}
	
	if ( (statDa <= valcritleft) || (statDa >= valcritright)) decision[i] = 1; else decision[i] = 0; 
	
      }
      
      if (usecrit[0] == 1) { // We use the provided critical values
	if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	} else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	} else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
      } else {
	//     if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }

// If applicable, we free the unused array of pointers
  delete[] xs;

}

// We return
  return;
   
        
  }
  
}

/*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
POUR CALCULER LES BONNES VALEURS CRITIQUES 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



myfunc <- function(M=10^7) {

require(PoweR)
level=0.05
vectn = 500
vectlois = 2
vectstats = 22
  vectn.len <- length(vectn)
  lois.len <- length(vectlois)
  stats.len <- length(vectstats)
  loinames <- paste(rep(" ",20*lois.len), sep = "", collapse = "")
  statnames <- paste(rep(" ",20*stats.len), sep = "", collapse = "")

    decision.len <- stats.len*vectn.len*lois.len
    decision <- rep(0,decision.len)

res <- rep(0,M)

for (i in 1:M) {

   res[i] <- .C("calcpuiss",as.integer(1),as.integer(vectlois),as.integer(lois.len),as.integer(vectn),as.integer(vectn.len),as.integer(vectstats),as.integer(stats.len),decision=as.integer(decision),as.integer(decision.len),level=as.double(level),loinames=as.character(loinames),statnames=as.character(statnames),PACKAGE="PoweR")$level

}

return(res)

}

require(parallel)

system.time({
nbclus <- 6

cl <- makeCluster(nbclus, type = "MPI") 
out <- clusterCall(cl, myfunc)

})

stopCluster(cl)



print(round(quantile(unlist(out),c(0.025,0.05,0.95,0.975)),3))



 */
