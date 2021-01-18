// g++ -c powcompeasy.cpp -o calcpuiss.o -I"/usr/lib/R/include"
// g++ -shared -o powcompeasy.so powcompeasy.o -I"/usr/lib/R/include" -L"/usr/lib" -lR


// WARNING: getname[0] should never be set to 1 in this file because it has the side effect to change the values of the parameters of either the laws or the stats.

#include <iostream>
// // using namespace std; // Pierre modified this on 22/12/2015

#include <R.h>
#include "Rmath.h"
#include <R_ext/Rdynload.h>

#include "laws-stats/def-laws-stats.cpp"

#include "models/def-models.cpp"

extern "C" {

static char *sfunction;

  // Generation of the sample
  
  int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale) {
    
    int i, n = xlen[0];
    double meanX = 0.0, varX = 0.0, sdX;

    (*lawfunc[law-1])(xlen, x, name, getname, params, nbparams, setseed);

    if (scale[0] == 1) {
      for (i = 0; i <= (n - 1); i++) meanX = meanX + x[i];
      meanX = meanX / (double)n;
      for (i = 0; i <= (n - 1); i++) varX = varX + R_pow(x[i], 2.0);
      varX = ((double)n) * (varX / (double)n - R_pow(meanX, 2.0)) / (double)(n - 1); 
      sdX = sqrt(varX);
      if (center[0] == 1) {
	for (i = 0; i <= (n - 1); i++) x[i] = (x[i] - meanX) / sdX;
      } else {
	for (i = 0; i <= (n - 1); i++) x[i] = x[i] / sdX;
      }
    } else {
      if (center[0] == 1) {
	for (i = 0; i <= (n - 1); i++) meanX = meanX + x[i];
	meanX = meanX / (double)n;
	for (i = 0; i <= (n - 1); i++) x[i] = x[i] - meanX;
      }
    }
    
    return(1);
  
  }
  
  // Computation of the test statistic
  
  void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {
    
    //    int i;

    //    double *xcopy; // POURQUOI J'AI FAIT CA??? IL FAUT LE FAIRE AUSSI POUR statistic, pvalue, etc ...?
    //    xcopy = new double[*(xlen+0)];
    //    for (i=1;i<=*(xlen+0);i++) *(xcopy+i-1) = *(x+i-1);    

    //    (*statfunc[stat-1])(xcopy,xlen,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparamstat);
    (*statfunc[stat-1])(x,xlen,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparamstat);
    
    //    delete[] xcopy;

    return;
    
  }

  // Models
  
  int model(int modelnum, char** funclist, double *thetavec, double *xvec, int *xlen, double *x, int *p, int *np) {
    
    if (modelnum == 0) { // Appel d'une fonction en R

      void call_R(char *func, long nargs, void **arguments, char **modes,long *lengths, char **names, long nres, char **results);

      int i;
      sfunction = funclist[0];
      
      long nargs=3; // Nombre d'arguments
      
      void **arguments;
      arguments = new void*[3];
      arguments[0] = (char*)&(x[0]);
      arguments[1] = (char*)&(thetavec[0]);
      arguments[2] = (char*)&(xvec[0]);
      
      char **modes; // Les modes de chaque argument
      modes = new char*[3];
      modes[0] = (char *)"double"; // Mode des éléments du 1er argument
      modes[1] = (char *)"double"; // Mode des éléments du 2ème argument
      modes[2] = (char *)"double"; // Mode des éléments du 3ème argument
      
      long *lengths; // Longueurs de chaque argument
      lengths = new long[3];
      lengths[0] = (long)(xlen[0]); // Longueur du 1er argument
      lengths[1] = (long)(p[0]); // Longueur du 2ème argument
      lengths[2] = (long)(np[0]); // Longueur du 3ème argument
      
      char **names; // Noms des arguments de la fonction provenant de funclist
      names = (char**)0;
      
      long nres=1; // Longeur de la liste de valeurs renvoyées par la fonction R dont le nom est stocké dans funclist
      
      char **results; // Les résultats
      results = new char*[1];
      results[0] = new char[xlen[0]];
      
      
      call_R(sfunction,nargs,arguments,modes,lengths,names,nres,results);
      
  
      for (i=0;i<=(xlen[0]-1);i++) {
	
		x[i] = ((double*)results[0])[i];
	
      }
      
      
      //On libere de la memoire
      delete[] arguments;
      delete[] modes;
      delete[] results;
      delete[] names;
      delete[] lengths;
      

    } else { // Appel d'une fonction en C

      (*modelfunc[modelnum-1])(xlen,x);
      
    }
    
    return(1);
  }

  
  // Computation of the power of the test statistic
  void powcompeasy(int *M, double *params, int *ncolparams, int *decision, int *decisionlen, //char **lawnames, char **statnames, 
		   int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale) {
    
    int gensample(int law, int *xlen, double *x, char **name1, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name2, int *getname, double *statistic, 
		     int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
    int model(int modelnum, char** funclist, double *thetavec, double *xvec, int *xlen, double *x, int *p, int *np);
    
    double *statistic, *pvalue; // POUR L'INSTANT JE N'EN FAIT RIEN DE statistic et de pvalue!! Si je veux les récupérer dans R il faudra faire des modifs!! A voir ...
    int *pvalcomp;
    statistic = new double[1];
    pvalue = new double[1];
    pvalcomp = new int[1];
    statistic[0] = 0.0;
    pvalue[0] = 0.0;
    pvalcomp[0] = 1;
    
    int i, row, n, law, stat, j, *xlen, *alter, *usecrit, *getname;
    double *x, *level, *critvalL, *critvalR, *parlaw, *parstat;
    char **name1, **name2;
    name1 = new char*[50];
    name2 = new char*[50];
    for (i = 0; i < 50; i++) {
      name1[i] =  new char[1];
      name2[i] =  new char[1];
      name1[i][0] = ' ';
      name2[i][0] = ' ';
    }
    getname = new int[1];
    int *decisiontmp, *nblevel, *nbparlaw, *nbparstat;
    decisiontmp = new int[1];
    decisiontmp[0] = 0;
    nblevel = new int[1];
    nblevel[0] = 1;
    
    int *setseed;
    setseed = new int[1];
    setseed[0] = 0;
    GetRNGstate();
    
    for (row = 1; row <= decisionlen[0]; row++) { // numéro de ligne de params considérée comme une matrice
      
      // lecture des valeurs des différentes colonnes de params sur la ligne row
      n = (int)params[(ncolparams[0]) * (row - 1) + 0];
      law = (int)params[(ncolparams[0]) * (row - 1) + 1];
      stat = (int)params[(ncolparams[0]) * (row - 1) + 2];
      level = new double[1];
      level[0] = params[(ncolparams[0]) * (row - 1) + 3];
      critvalL = new double[1];
      critvalL[0] = params[(ncolparams[0]) * (row - 1) + 4];
      critvalR = new double[1];
      critvalR[0] = params[(ncolparams[0]) * (row - 1) + 5];
      alter = new int[1];
      alter[0] = (int)params[(ncolparams[0]) * (row - 1) + 6];
      usecrit = new int[1];
      usecrit[0] = (int)params[(ncolparams[0]) * (row - 1) + 7];
      nbparlaw = new int[1];
      nbparlaw[0] = (int)params[(ncolparams[0]) * (row - 1) + 8];
      parlaw = new double[4];   // j'ai considéré que 4 paramètres étaient suffisants pour chaque loi. Si jamais on propose une loi à 5 paramètres, il faudra modifier le code
      parlaw[0] = params[(ncolparams[0]) * (row - 1) + 9];
      parlaw[1] = params[(ncolparams[0]) * (row - 1) + 10];
      parlaw[2] = params[(ncolparams[0]) * (row - 1) + 11];
      parlaw[3] = params[(ncolparams[0]) * (row - 1) + 12];
      nbparstat = new int[1];
      nbparstat[0] = (int)params[(ncolparams[0]) * (row - 1) + 13];
      parstat = new double[nbparstat[0]];  
      for (i = 0; i < nbparstat[0]; i++) parstat[i] = params[(ncolparams[0]) * (row - 1) + 14 + i];

      x = new double[n];
      for (j = 0; j < n; j++) x[j] = 0.0;
      xlen = new int[1];
      xlen[0] = n;

      getname[0] = 0;
      if (M[0] > 1) {
	for (i = 1; i <= M[0]; i++) { // on part la simul, sans refaire la 1ere iteration!
	  
	  gensample(law, xlen, x, name1, getname, parlaw, nbparlaw, setseed, center, scale); // on génère l'échantillon
	  model(modelnum[0], funclist, thetavec, xvec, xlen, x, p, np);  // on applique le modèle
	  
	  statcompute(stat, x, xlen, level, nblevel, name2, getname, statistic, pvalcomp, pvalue, critvalL, critvalR, usecrit, alter, decisiontmp, parstat, nbparstat);
	  decision[row-1] = decision[row-1] + decisiontmp[0];
	  
	}
      }

     // We retrieve the default values of parameter laws and parameter stats used 
      params[(ncolparams[0]) * (row - 1) + 8]  = (int)nbparlaw[0];
      params[(ncolparams[0]) * (row - 1) + 9]  = parlaw[0];
      params[(ncolparams[0]) * (row - 1) + 10] = parlaw[1];
      params[(ncolparams[0]) * (row - 1) + 11] = parlaw[2];
      params[(ncolparams[0]) * (row - 1) + 12] = parlaw[3];
      params[(ncolparams[0]) * (row - 1) + 13] = (int)nbparstat[0];
      for (i = 0; i < nbparstat[0]; i++) params[(ncolparams[0]) * (row - 1) + 14 + i] = parstat[i];
    
      //On libere de la memoire
      delete[] x;
      delete[] xlen;
      delete[] level;
      delete[] critvalL;
      delete[] critvalR;
      delete[] alter;
      delete[] usecrit;
      delete[] nbparlaw;
      delete[] parlaw;
      delete[] nbparstat;
      delete[] parstat;
      
    }
    

    //On libere de la memoire
    for (i = 0; i < 50; i++) {
      delete[] name1[i];
      delete[] name2[i];
    }
    delete[] name1;
    delete[] name2;
    delete[] statistic;
    delete[] pvalue;
    delete[] pvalcomp;
    delete[] decisiontmp;
    delete[] nblevel;
    delete[] setseed;
    delete[] getname;    

    PutRNGstate();
    return;
    
  }

  // computation of all the statistic values in order to obtain the critical values
  void compquantc(int *n, int *law, int *stat, int *M, double *statvec, //char **lawname, char **statname, 
		 int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale) {

    int gensample(int law, int *xlen, double *x, char **name1, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp,double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
    int model(int modelnum, char** funclist, double *thetavec, double *xvec, int *xlen, double *x, int *p, int *np);
    
    double *statistic, *pvalue, *level, *critvalL, *critvalR; // POUR L'INSTANT JE N'EN FAIT RIEN DE statistic et de pvalue!! Si je veux les récupérer dans R il faudra faire des modifs!! A voir ...
    int *pvalcomp;
    level = new double[1];
    statistic = new double[1];
    pvalue = new double[1];
    pvalcomp = new int[1];
    critvalL = new double[1];
    critvalR = new double[1];
    *(level+0) = 0.0;
    *(statistic+0) = 0.0;

	// We set pvalcomp[0] = 0 so the function compquant() doesn't compute p-value each time it is called
	// if you want to retrieve these p-values, please change pvalue[0] to 0
    pvalue[0] = 0.0;
    pvalcomp[0] = 0;
	
	// end
	
    *(critvalL+0) = 0.0;
    *(critvalR+0) = 0.0;

    int *usecrit, *alter;
    usecrit = new int[1];
    alter = new int[1];
    *(usecrit+0) = 0;
    *(alter+0) = 0;

    int *decisiontmp, *nblevel;
    decisiontmp = new int[1];
    decisiontmp[0] = 0;
    nblevel = new int[1];
    nblevel[0] = 1;

    int i, j, *getname;
    getname = new int[1];

    double *x;

    char **name1, **name2;
    name1 = new char*[50];
    name2 = new char*[50];
    for (i = 0; i < 50; i++) {
      name1[i] =  new char[1];
      name2[i] =  new char[1];
      name1[i][0] = ' ';
      name2[i][0] = ' ';
    }

    int *setseed;
    setseed = new int[1];
    setseed[0] = 0;
    GetRNGstate();
	      
    getname[0] = 0;	
    for (i=1;i<=*(M+0);i++) {
      
      x = new double[*(n+0)];
      for (j=1;j<=*(n+0);j++) *(x+j-1) = 0.0;
      gensample(*(law+0),n,x,name1,getname,parlaw,nbparlaw,setseed,center,scale);
      model(modelnum[0],funclist,thetavec,xvec,n,x,p,np);   // on applique le modèle
    
      statcompute(*(stat+0), x, n, level, nblevel, name2, getname, statistic, pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decisiontmp,parstat,nbparstat);
    
      *(statvec+i-1) = *(statistic+0);      

      //On libere de la memoire
      delete[] x;
      	      
      
	      }

    //On libere de la memoire
    for (i = 0; i < 50; i++) {
      delete[] name1[i];
      delete[] name2[i];
    }
    delete[] name1;
    delete[] name2;
    delete[] level;
    delete[] statistic;
    delete[] pvalue;
    delete[] pvalcomp;
    delete[] critvalL;
    delete[] critvalR;
    delete[] usecrit;
    delete[] alter;
    delete[] decisiontmp;
    delete[] nblevel;
    delete[] setseed;
    delete[] getname;    

    PutRNGstate();
    return;
    
  }


  // Computation of the power of the test statistic
  void powcompfast(int *M, int *vectlaws, int *lawslen, int *vectn, int *vectnlen, int *vectstats, int *statslen, int *decision, int *decisionlen, double *level, int *nblevel,
		   double *critvalL, double *critvalR, int *usecrit, int *alter, int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, 
		   double *thetavec, double *xvec, int *p, int *np, int *center, int *scale, int *compquant) {


    // Warning: When compquant[0] == 1, critvalL should be initialized with lawslen[0] * M[0] * vectnlen[0] * statslen[0] double values since it will contain (when output)
    // all the test statistics generated. Most probably, lawslen[0] should be equal to 1.

    // Il faut permuter M et law ci-dessous et aussi dans les fichiers R concernés (ainsi que dans la fonction print)!!

    int gensample(int law, int *xlen, double *x, char **name1, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name2, int *getname, double *statistic, 
		     int *pvalcomp,double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
    int model(int modelnum, char** funclist, double *thetavec, double *xvec, int *xlen, double *x, int *p, int *np);
   
    double *statistic, *pvalue; // POUR L'INSTANT JE N'EN FAIT RIEN DE pvalue!! Si je veux les récupérer dans R il faudra faire des modifs!! A voir ...
    int *pvalcomp;
    statistic = new double[1];
    pvalue = new double[1];
    pvalcomp = new int[1];
    statistic[0] = 0.0;
    pvalue[0] = 0.0;
    if (compquant[0] == 1) pvalcomp[0] = 0; else pvalcomp[0] = 1;

    int indtmp = 0;

    int i, n, law, stat, j, *xlen, *getname;
    double *x;
    getname = new int[1];

    char **name1, **name2;
    name1 = new char*[50];
    name2 = new char*[50];
    for (i = 0; i < 50; i++) {
      name1[i] =  new char[1];
      name2[i] =  new char[1];
      name1[i][0] = ' ';
      name2[i][0] = ' ';
    }

    int *decisiontmp;
    decisiontmp = new int[nblevel[0]];
    for (i = 0; i < nblevel[0]; i++) decisiontmp[i] = 0;

    int *setseed;
    setseed = new int[1];
    setseed[0] = 0;
    GetRNGstate();


    int maxn = vectn[0]; // maximum of vectn values
    for (i = 1; i < vectnlen[0]; i++) if (maxn < vectn[i]) maxn = vectn[i];

    x = new double[maxn]; // on reserve un vecteur de taille (max des n dans vectn)
    for (j = 0; j < maxn; j++) x[j] = 0.0;
    xlen = new int[1];

    int *altertmp, *usecrittmp, *nbparlawtmp, *nbparstattmp;
    altertmp = new int[1];
    usecrittmp = new int[1];
    nbparlawtmp = new int[1];
    nbparstattmp = new int[1];
	
    double *critvalLtmp, *critvalRtmp, *parlawtmp, *parstattmp;
    critvalLtmp = new double[nblevel[0]];
    critvalRtmp = new double[nblevel[0]];
    parlawtmp = new double[4];
	
    int t;
    //    int m;
    int kmax = 0; 		// kmax = parstats.len.max in powcomp-fast.R
    for (t = 0; t < statslen[0]; t++) {
      if (kmax <= nbparstat[t]) {
	kmax = nbparstat[t];
      } //else kmax = kmax;
    }
    if (kmax == 0) kmax = 1;
    parstattmp = new double[kmax];

    int stlen1, stlen2;

    getname[0] = 0;	

    for (law = 0; law < lawslen[0]; law++) {

      for (i = 1; i <= M[0]; i++) {

	xlen[0] = maxn; // on génère un échantillon de taille maxn de loi law 

	nbparlawtmp[0] = nbparlaw[law];
	parlawtmp[0] = parlaw[0 + 4 * law];
	parlawtmp[1] = parlaw[1 + 4 * law];
	parlawtmp[2] = parlaw[2 + 4 * law];
	parlawtmp[3] = parlaw[3 + 4 * law];
	
	gensample(vectlaws[law], xlen, x, name1, getname, parlawtmp, nbparlawtmp, setseed, center, scale);
	model(modelnum[0], funclist, thetavec, xvec, xlen, x, p, np);   // on applique le modèle
	
	if (i == 1) {
	  nbparlaw[law] = nbparlawtmp[0];
	  parlaw[0 + 4 * law] = parlawtmp[0];
	  parlaw[1 + 4 * law] = parlawtmp[1];
	  parlaw[2 + 4 * law] = parlawtmp[2];
	  parlaw[3 + 4 * law] = parlawtmp[3];
	}

	for (n = 0; n < vectnlen[0]; n++) {
	  
	  xlen[0] = vectn[n]; // permet de ne prendre que la portion (au début) de x de taille vectn[n-1]
	  
	  stlen1 = 0; stlen2 = 0;
	  for (stat = 0; stat < statslen[0]; stat++) {
	    
	    altertmp[0] = alter[stat];
	    
	    nbparstattmp[0] = nbparstat[stat];
	    for (t = 0; t < nbparstattmp[0]; t++) {
	      parstattmp[t] = parstat[t + stlen1];
	    }
	    stlen1 = stlen1 + nbparstattmp[0];
	
	    usecrittmp[0] = usecrit[n + vectnlen[0] * stat];
	    for (j = 0; j < nblevel[0]; j++) {
	      indtmp = n + nblevel[0] * vectnlen[0] * stat + vectnlen[0] * j;
	      critvalLtmp[j] = critvalL[indtmp];
	      critvalRtmp[j] = critvalR[indtmp];	  
	    }

	    // decisiontmp est de longueur nblevel[0]
	    statcompute(vectstats[stat], x, xlen, level, nblevel, name2, getname, statistic, pvalcomp, pvalue, critvalLtmp, critvalRtmp, usecrittmp, altertmp, decisiontmp, parstattmp, nbparstattmp);
	    if (compquant[0] == 1) {
	      critvalL[i + M[0] * n + M[0] * vectnlen[0] * stat - 1] = statistic[0];
	    }


	    if (i == 1 && law == 0 && n == 0) {
	      nbparstat[stat] = nbparstattmp[0];
	      for (t = 0; t < nbparstattmp[0]; t++) {
		parstat[t + stlen2] = parstattmp[t];
	      }
	      stlen2 = stlen2 + nbparstattmp[0];
	    }

	    if (compquant[0] == 0) {
	      for (j = 0; j < nblevel[0]; j++) {
		indtmp = stat + statslen[0] * n + statslen[0] * vectnlen[0] * law + statslen[0] * vectnlen[0] * lawslen[0] * j;
		decision[indtmp] = decision[indtmp] + decisiontmp[j];
	      }
	    }
	    
	  }
	  
	}
	
      }
      
    }

    //On libere de la memoire
    for (i = 0; i < 50; i++) {
      delete[] name1[i];
      delete[] name2[i];
    }
    delete[] x;
    delete[] name1;
    delete[] name2;
    delete[] statistic;
    delete[] pvalue;
    delete[] pvalcomp;
    delete[] decisiontmp;
    delete[] xlen;
    delete[] altertmp;
    delete[] usecrittmp;
    delete[] critvalLtmp;
    delete[] critvalRtmp;
    delete[] setseed;
    delete[] nbparlawtmp;
    delete[] parlawtmp;
    delete[] nbparstattmp;
    delete[] parstattmp;
    delete[] getname;

    PutRNGstate();
    return;

  }
  
  
  // Function calcFx() version C++
  void calcfx(double *pvals, int *pvalslen, double *xi, int *xilen, double *fx) {
  
	int i, j;
	int tmp;
	
	for (j=0; j<xilen[0]; j++) {
	
	  tmp = 0;
	
	  for (i=0; i<pvalslen[0]; i++) {
	  
	    if (pvals[i] <= xi[j]) tmp = tmp + 1;
	   
	  }
	
	  fx[j] = (double)tmp/(double)pvalslen[0];
	
	}
	
	return;
  
  }









  
   
  // Computation of the p-values matrix
  void matrixpval(int* N, int* lawindex, int *xlen, int *nbparams, double *parlaw, int *statindices, int* nbstats, int *altervec, double *parstatmultvec, int *nbparstatvec, double *res, int *center, int *scale) {

    int i, j, k;

    int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);

    double *x; 
    x = new double[*(xlen+0)];
    for (i=0;i<xlen[0];i++) x[i] = 0.0; 

    double *params;
    params = new double[4];
    if (nbparams[0] == 0) {
      params[0] = 0.0;
      params[1] = 0.0;
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 1) {
      params[0] = parlaw[0];
      params[1] = 0.0;
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 2) {
      params[0] = parlaw[0];
      params[1] = parlaw[1];
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 3) {
      params[0] = parlaw[0];
      params[1] = parlaw[1];
      params[2] = parlaw[2];
      params[3] = 0.0;
    }

    char **name;          // name = rep(" ", 50)
    name = new char*[50];
    for (i = 0; i < 50; i++) {
      name[i] =  new char[1];
      name[i][0] = ' ';
    }

    int *getname;
    getname = new int[1];
    getname[0] = 0;

    int *setseed;
    setseed = new int[1];
    setseed[0] = 1;
    GetRNGstate();

    int statindex, *nblevel, *usecrit, *decision;
    nblevel = new int[1];
    usecrit = new int[1];
    decision = new int[1];
    double *level, *critvalL, *critvalR, *statistic, *pvalue;
    int *pvalcomp;
    level = new double[1];
    critvalL = new double[1];
    critvalR = new double[1];
    statistic = new double[1];
    pvalue = new double[1];
    pvalcomp = new int[1];

    int *alter;
    alter= new int[1];

    int * nbparamstat;
    nbparamstat = new int[1];


    int cmpt;

    level[0] = 0.05;
    nblevel[0] = 1;
    usecrit[0] = 0;
    critvalL[0] = 0.0;
    critvalR[0] = 0.0;
    statistic[0] = 0.0;
    pvalue[0] = 0.0;
    decision[0] = 0;

    for (j=0;j<N[0];j++) {
      
      gensample(lawindex[0],xlen,x,name,getname,params,nbparams,setseed,center,scale);  // We retrieve x
      
      cmpt = 0;

      for (i=0;i<nbstats[0];i++) {

	pvalcomp[0] = 1;

	statindex = statindices[i];	    

	alter[0] = altervec[i];

	nbparamstat[0] = nbparstatvec[i];

	if (nbparamstat[0]>0) {

	  double * paramstat;
	  paramstat = new double[nbparamstat[0]];    
	  for (k=0;k<nbparamstat[0];k++) paramstat[k] = parstatmultvec[cmpt+k];
	  cmpt = cmpt + nbparamstat[0];
	  (*statfunc[statindex-1])(x,xlen,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparamstat);
	  delete[] paramstat;
	} else {
	  (*statfunc[statindex-1])(x,xlen,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,(double*)0,nbparamstat);
	}

	res[i*N[0] + j]  = pvalue[0];
	// Avant je pensais pouvoir retourner des NA là dedans et m'en servir ?? 
	// Devrais-je retourner pvalue[0] = 2 si jamais pvalcomp[0] = 0 ??
	//	if (pvalcomp[0] == 1) res[i*N[0] + j]  = pvalue[0]; else res[i*N[0] + j]  = 2;

      }
    }

    delete[] x;
    delete[] params;
    for (i=0;i<50;i++) {
      delete[] name[i];
    }
    delete[] name;
    delete[] getname;
    delete[] setseed;
    delete[] nblevel;
    delete[] usecrit;
    delete[] decision;
    delete[] level;
    delete[] critvalL;
    delete[] critvalR;
    delete[] statistic;
    delete[] pvalue;
    delete[] pvalcomp;
    delete[] alter;
    delete[] nbparamstat;

    PutRNGstate();
    return;

  }
  


  // Computation of the p-values matrix using Monte-Carlo
  void matrixpvalMC(int *n, int *lawindex, int* nbstats, int *M, int *statindices, int *nbparstatvec, double *parstatmultvec, char** funclist, int *N, int *nulldist, int *nbparams, int *altervec, double *parstat, int *nbparstat, double *res, int *center, int *scale) {

    void compquantc(int *n, int *law, int *stat, int *M, double *statvec,// char **lawname, char **statname, 
int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale);
    int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);


    int i, j, k;

    int *stat;
    stat = new int[1];

    double *statvec;
    statvec = new double[M[0]];
    for (i=0;i<M[0];i++) statvec[i] = 0.0;



    char **lawname;         // lawname=paste(rep(" ",50), sep = "", collapse = "")
    lawname = new char*[1];
    lawname[0] = new char[50];


    char **statname;          // statname = paste(rep(" ",50), sep = "", collapse = "")
    statname = new char*[1];
    statname[0] = new char[50];

    int *nbparlaw;
    nbparlaw = new int[1];
    nbparlaw[0] = 0;

    double *parlaw;
    parlaw = new double[4];
    parlaw[0] = 0.0; parlaw[1] = 0.0; parlaw[2] = 0.0; parlaw[3] = 0.0; 


    int *nbparamstat;
    nbparamstat = new int[1];

    int cmpt;

    int *modelnum;
    modelnum = new int[1];
    modelnum[0] = 1;

    double *thetavec;
    thetavec = new double[1];
    thetavec[0] = 0.0;

    double *xvec;
    xvec = new double[1];
    xvec[0] = 0.0;

    int *p;
    p = new int[1];
    p[0] = 1;

    int *np;
    np = new int[1];
    np[0] = 1;


    double *liststat;
    liststat = new double[M[0]*nbstats[0]];

      cmpt = 0;

    for (i=0;i<nbstats[0];i++) {

      stat[0] = statindices[i];

      nbparamstat[0] = nbparstatvec[i];
      
      if (nbparamstat[0]>0) {	
	double * paramstat;
	paramstat = new double[nbparamstat[0]];    
	for (k=0;k<nbparamstat[0];k++) paramstat[k] = parstatmultvec[cmpt+k];
	cmpt = cmpt + nbparamstat[0];
	compquantc(n,lawindex,stat,M,statvec,//lawname,statname,
		  nbparlaw,parlaw,nbparamstat,paramstat,modelnum,funclist,thetavec,xvec,p,np,center,scale);
	delete[] paramstat;
      } else {
	compquantc(n,lawindex,stat,M,statvec,// lawname,statname,
		  nbparlaw,parlaw,nbparamstat,(double*)0,modelnum,funclist,thetavec,xvec,p,np,center,scale);
      }

      for (j=0;j<M[0];j++) liststat[i*M[0] + j] = statvec[j];

    }



    double *x; 
    x = new double[n[0]];
    for (i=0;i<n[0];i++) x[i] = 0.0; 


    double *params;
    params = new double[4];
    if (nbparams[0] == 0) {
      params[0] = 0.0;
      params[1] = 0.0;
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 1) {
      params[0] = parlaw[0];
      params[1] = 0.0;
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 2) {
      params[0] = parlaw[0];
      params[1] = parlaw[1];
      params[2] = 0.0;
      params[3] = 0.0;
    }
    if (nbparams[0] == 3) {
      params[0] = parlaw[0];
      params[1] = parlaw[1];
      params[2] = parlaw[2];
      params[3] = 0.0;
    }


    char **name;          // name = rep(" ", 50)
    name = new char*[50];
    for (i = 0; i < 50; i++) {
      name[i] =  new char[1];
      name[i][0] = ' ';
    }

    int *getname;
    getname = new int[1];
    getname[0] = 0;

    int *setseed;
    setseed = new int[1];
    setseed[0] = 1;
    GetRNGstate();


    int *nblevel, *usecrit, *decision;
    nblevel = new int[1];
    usecrit = new int[1];
    decision = new int[1];
    double *level, *critvalL, *critvalR, *statistic, *pvalue;
    int *pvalcomp;
    level = new double[1];
    critvalL = new double[1];
    critvalR = new double[1];
    statistic = new double[1];
    pvalue = new double[1];
    pvalcomp = new int[1];




    double *liststati;
    liststati = new double[M[0]];

    double q2;

    double meanstatsup, meanstatinf;

    int *alter;
    alter = new int[1];

    for (j=0;j<N[0];j++) {

      gensample(nulldist[0],n,x,name,getname,params,nbparams,setseed,center,scale);  // We retrieve x

      cmpt = 0;
      for (i=0;i<nbstats[0];i++) {

	level[0] = 0.05;
	nblevel[0] = 1;
	usecrit[0] = 0;
	critvalL[0] = 0.0;
	critvalR[0] = 0.0;
	statistic[0] = 0.0;
	pvalue[0] = 0.0;
	pvalcomp[0] = 0;	// We set pvalcomp[0] = 0 so the function statxxx doesn't compute p-value each time it is called

	decision[0] = 0;
	alter[0] = altervec[i];

	if (nbparamstat[0]>0) {	
	  double * paramstat;
	  paramstat = new double[nbparamstat[0]];    
	  for (k=0;k<nbparamstat[0];k++) paramstat[k] = parstatmultvec[cmpt+k];
	  cmpt = cmpt + nbparamstat[0];
	  (*statfunc[statindices[i]-1])(x,n,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparstat);    // We obtain one test statistic
	  delete[] paramstat;
	} else {
	  (*statfunc[statindices[i]-1])(x,n,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,(double*)0,nbparstat);    // We obtain one test statistic
	}		
	
	for (k=0;k<M[0];k++) liststati[k] = liststat[i*M[0] + k];	
	R_rsort(liststati,M[0]); 
	if(M[0] % 2 == 0) {		// check if M is divisible by 2
	  q2 = (liststati[M[0]/2-1] + liststati[M[0]/2])/2.0;
	} else {
	  q2 = liststati[M[0]/2];
	}
			
	// calcul pvalueMC by method of Fisher
	meanstatsup = 0.0;
	meanstatinf = 0.0;
	for (k=0;k<M[0];k++) {
	  if (liststati[k] >= statistic[0]) meanstatsup = meanstatsup + 1;
	  if (liststati[k] <= statistic[0]) meanstatinf = meanstatinf + 1;
	}
	meanstatsup = meanstatsup/M[0];
	meanstatinf = meanstatinf/M[0];			
	if (alter[0] == 0) {	
	  if (statistic[0] >= q2) {
	    pvalue[0] = 2*meanstatsup;
	  } else {
	    pvalue[0] = 2*meanstatinf;
	  }
	}		
	if ((alter[0] == 1) | (alter[0] == 4)) {
	  pvalue[0] = meanstatinf;
	}	
	if ((alter[0] == 2) | (alter[0] == 3)) {
	  pvalue[0] = meanstatsup;
	}
	      	
	res[i*N[0] + j] = pvalue[0];
	
      }

    }

    delete[] stat;
    delete[] statvec;
    delete[] lawname[0];
    delete[] lawname;
    delete[] statname[0];
    delete[] statname;
    delete[] nbparlaw;
    delete[] parlaw;
    delete[] nbparamstat;
    delete[] modelnum;
    delete[] thetavec;
    delete[] xvec;
    delete[] p;
    delete[] np;
    delete[] x;
    delete[] getname;
    delete[] params;
    delete[] setseed;
    for (i=1;i<=50;i++) {
      delete[] *(name+i-1);
    }
    delete[] name;
    delete[] nblevel;
    delete[] usecrit;
    delete[] decision;
    delete[] level;
    delete[] critvalL;
    delete[] critvalR;
    delete[] statistic;
    delete[] pvalue;
    delete[] pvalcomp;
    delete[] liststati;
    delete[] liststat;
    delete[] alter;

    PutRNGstate();
    return;

  }

#include "Rcpp/calcpuissRcpp.cpp"

  
#include "laws-stats/register.cpp"


}


