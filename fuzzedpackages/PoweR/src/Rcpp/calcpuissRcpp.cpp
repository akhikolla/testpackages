#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;


// gensample <- function(law,n,params) {return(.Call("gensampleRcpp",law,n,params,as.character(match.call()[2]),PACKAGE="PoweR"))}

extern "C" {


  SEXP gensampleRcpp2( Function rlawfunc, IntegerVector n, NumericVector params , int nbparams, const std::string lawname, IntegerVector center, IntegerVector scale) {
  //	 Rcpp::RNGScope __rngScope;
  if (nbparams > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase this value.");
  Rcpp::List result;
  if (nbparams == 0) {
    result = Rcpp::List::create(Rcpp::Named("sample") = rlawfunc(n),
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = R_NilValue);
  } else if (nbparams == 1) {
    result = Rcpp::List::create(Rcpp::Named("sample") = rlawfunc(n,params[0]),
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = params);

  } else if (nbparams == 2) {
    result = Rcpp::List::create(Rcpp::Named("sample") = rlawfunc(n,params[0],params[1]),
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = params);

  } else if (nbparams == 3) {
    result = Rcpp::List::create(Rcpp::Named("sample") = rlawfunc(n,params[0],params[1],params[2]),
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = params);

  } else if (nbparams == 4) {
    result = Rcpp::List::create(Rcpp::Named("sample") = rlawfunc(n,params[0],params[1],params[2],params[3]),
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = params);
  } else {
    result = Rcpp::List::create(Rcpp::Named("sample") = 0,
			      Rcpp::Named("law.name") = lawname,
			      Rcpp::Named("law.pars") = params);
  }

  int i;
  double meanX = 0.0, varX = 0.0, sdX;

  NumericVector x = as<NumericVector>(result["sample"]);

    if (scale[0] == 1) {
      for (i = 0; i <= (n[0] - 1); i++) meanX = meanX + x[i];
      meanX = meanX / (double)n[0];
      for (i = 0; i <= (n[0] - 1); i++) varX = varX + R_pow(x[i], 2.0);
      varX = ((double)n[0]) * (varX / (double)n[0] - R_pow(meanX, 2.0)) / (double)(n[0] - 1); 
      sdX = sqrt(varX);
      if (center[0] == 1) {
	for (i = 0; i <= (n[0] - 1); i++) x[i] = (x[i] - meanX) / sdX;
      } else {
	for (i = 0; i <= (n[0] - 1); i++) x[i] = x[i] / sdX;
      }
    } else {
      if (center[0] == 1) {
	for (i = 0; i <= (n[0] - 1); i++) meanX = meanX + x[i];
	meanX = meanX / (double)n[0];
	for (i = 0; i <= (n[0] - 1); i++) x[i] = x[i] - meanX;
      }
    }

    result["sample"] = x;

  return result;
	
}

RcppExport SEXP gensampleRcpp(SEXP rlawfuncSEXP, SEXP nSEXP, SEXP paramsSEXP , SEXP nbparamsSEXP, SEXP lawnameSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
  BEGIN_RCPP
    //   Rcpp::RNGScope __rngScope; // useful? takes time ... and is used for setting seeds ..
  Function rlawfunc = Rcpp::as<Function >(rlawfuncSEXP);
  IntegerVector n = Rcpp::as<IntegerVector >(nSEXP);
  NumericVector params = Rcpp::as<NumericVector >(paramsSEXP);
  int nbparams = Rcpp::as<int >(nbparamsSEXP);
  std::string lawname = Rcpp::as<std::string >(lawnameSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  SEXP __result = gensampleRcpp2(rlawfunc, n, params, nbparams, lawname, center, scale);
  return Rcpp::wrap(__result);
  END_RCPP
    }


  SEXP statcomputeRcpp2( Function rstatfunc, SEXP ech, SEXP levels, SEXP usecrit, SEXP critvalL, SEXP critvalR) {
    //	 Rcpp::RNGScope __rngScope;
    List res = rstatfunc(ech,levels,usecrit,critvalL,critvalR);
	 return Rcpp::List::create(Rcpp::Named("statistic") = res["statistic"],
				   Rcpp::Named("pvalue") = res["pvalue"],
				   Rcpp::Named("decision") = res["decision"], // decision est de longueur nblevel[0]
				   Rcpp::Named("alter") = res["alter"],
				   Rcpp::Named("stat.pars") = res["stat.pars"],
				   Rcpp::Named("pvalcomp") = res["pvalcomp"],
				   Rcpp::Named("nbparstat") = res["nbparstat"]);
	
}



  RcppExport SEXP statcomputeRcpp(SEXP rstatfuncSEXP, SEXP echSEXP, SEXP levelsSEXP, SEXP usecritSEXP, SEXP critvalLSEXP, SEXP critvalRSEXP) {
  BEGIN_RCPP
    //   Rcpp::RNGScope __rngScope; // useful? takes time ... and is used for setting seeds ..
  Function rstatfunc = Rcpp::as<Function >(rstatfuncSEXP);
  NumericVector ech = Rcpp::as<NumericVector >(echSEXP);
  NumericVector levels = Rcpp::as<NumericVector >(levelsSEXP);
  IntegerVector usecrit = Rcpp::as<IntegerVector >(usecritSEXP);
  NumericVector critvalL = Rcpp::as<NumericVector >(critvalLSEXP);
  NumericVector critvalR = Rcpp::as<NumericVector >(critvalRSEXP);
  SEXP __result = statcomputeRcpp2(rstatfunc, ech, levels, usecrit, critvalL, critvalR);
  return Rcpp::wrap(__result);
  END_RCPP
    }


SEXP powcompeasyRcpp2(IntegerVector M, NumericVector params, IntegerVector ncolparams, IntegerVector decision, IntegerVector decisionlen,  
		      IntegerVector modelnum, List funclist, NumericVector thetavec, NumericVector xvec, IntegerVector p, IntegerVector np, List Rlaws, List Rstats, IntegerVector center, IntegerVector scale) {

  int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
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
    
    int i, k, row, n, law, stat, j, *xlen, *alter, *usecrit, *getname, *ptmp, *nptmp, *centertmp, *scaletmp;
    ptmp = new int[1];
    ptmp[0] = p[0];
    nptmp = new int[1];
    nptmp[0] = np[0];
    centertmp = new int[1];
    centertmp[0] = center[0];
    scaletmp = new int[1];
    scaletmp[0] = scale[0];
    double *x, *level, *critvalL, *critvalR, *parlaw, *parstat, *thetavectmp, *xvectmp;
    thetavectmp = new double[p[0]];
    for (i=1;i<=p[0];i++) thetavectmp[i-1] = 0.0;
    xvectmp = new double[np[0]];
    for (i=1;i<=np[0];i++) xvectmp[i-1] = 0.0;
    char **name1, **name2, **funclisttmp;
    funclisttmp = new char*[1];
    funclisttmp[0] = (char *)R_ExternalPtrAddr(funclist[0]);
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
    GetRNGstate(); Rcpp::RNGScope __rngScope;
    
    for (row=1;row<=decisionlen[0];row++) { // numéro de ligne de params considérée comme une matrice

      // lecture des valeurs des différentes colonnes de params sur la ligne row
      n = (int)params[(ncolparams[0])*(row-1)+0];
      law = (int)params[(ncolparams[0])*(row-1)+1];
      stat = (int)params[(ncolparams[0])*(row-1)+2];
      level = new double[1];
      level[0] = params[(ncolparams[0])*(row-1)+3];
      critvalL = new double[1];
      critvalL[0] = params[(ncolparams[0])*(row-1)+4];
      critvalR = new double[1];
      critvalR[0] = params[(ncolparams[0])*(row-1)+5];
      alter = new int[1];
      alter[0] = (int)params[(ncolparams[0])*(row-1)+6];
      usecrit = new int[1];
      usecrit[0] = (int)params[(ncolparams[0])*(row-1)+7];
      nbparlaw = new int[1];
      nbparlaw[0] = (int)params[(ncolparams[0])*(row-1)+8];
      parlaw = new double[4];   // j'ai considéré que 4 paramètres étaient suffisants pour chaque loi. Si jamais on propose une loi à 5 paramètres, il faudra modifier le code
      parlaw[0] = params[(ncolparams[0])*(row-1)+9];
      parlaw[1] = params[(ncolparams[0])*(row-1)+10];
      parlaw[2] = params[(ncolparams[0])*(row-1)+11];
      parlaw[3] = params[(ncolparams[0])*(row-1)+12];
      nbparstat = new int[1];
      nbparstat[0] = (int)params[(ncolparams[0])*(row-1)+13];
      parstat = new double[nbparstat[0]];  
      for (i=0;i<nbparstat[0];i++) parstat[i] = params[(ncolparams[0])*(row-1)+ 14 + i];

      x = new double[n];
      for (j=1;j<=n;j++) x[j-1] = 0.0;
      xlen = new int[1];
      xlen[0] = n;

      Function rlawfunc = Rcpp::as<Function >(Rlaws[row-1]);
      IntegerVector nSEXP(1);
      nSEXP[0] = n; 

      getname[0] = 0;
      if (law == 0) {
	if (M[0]>1) {
	  for (i=1;i<=M[0];i++) { // on part la simul, sans refaire la 1ere iteration!

	    List resultsample = gensampleRcpp2(rlawfunc, nSEXP, params, nbparlaw[0], "", center, scale);
	    NumericVector mysample = resultsample["sample"];
	    for (j=1;j<=n;j++) x[j-1] = mysample[j-1];
	    model(modelnum[0],funclisttmp,thetavectmp,xvectmp,xlen,x,ptmp,nptmp);  // on applique le modèle
	    
	    if (stat == 0) {
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[row-1],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      NumericVector mydecision = resultstat["decision"];
	      decision[row-1] = decision[row-1] + mydecision[0];
	    } else {
	      statcompute(stat, x, xlen, level, nblevel, name2, getname, statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decisiontmp,parstat,nbparstat);
	      decision[row-1] = decision[row-1] + decisiontmp[0];
	    }
	    
	  }
	}	
      } else {
	if (M[0]>1) {
	  for (i=1;i<=M[0];i++) { // on part la simul, sans refaire la 1ere iteration!
	    
	    gensample(law,xlen,x,name1,getname,parlaw,nbparlaw,setseed,centertmp,scaletmp); // on génère l'échantillon
	    model(modelnum[0],funclisttmp,thetavectmp,xvectmp,xlen,x,ptmp,nptmp);  // on applique le modèle
	    
	    if (stat == 0) {
	      NumericVector mysample (xlen[0]);
	      for (k=1;k<=M[0];k++) mysample[k-1] = x[k-1];
	      NumericVector levelSEXP = level[0];	      
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[row-1],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      NumericVector mydecision = resultstat["decision"];
	      decision[row-1] = decision[row-1] + mydecision[0];
	    } else {
	      statcompute(stat, x, xlen, level, nblevel, name2, getname, statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decisiontmp,parstat,nbparstat);
	      decision[row-1] = decision[row-1] + decisiontmp[0];
	    }
	    
	  }
	}    PutRNGstate();
      }

     // We retrieve the default values of parameter laws and parameter stats used 
      params[(ncolparams[0])*(row-1)+8]  = (int)nbparlaw[0];
      params[(ncolparams[0])*(row-1)+9]  = parlaw[0];
      params[(ncolparams[0])*(row-1)+10] = parlaw[1];
      params[(ncolparams[0])*(row-1)+11] = parlaw[2];
      params[(ncolparams[0])*(row-1)+12] = parlaw[3];
      params[(ncolparams[0])*(row-1)+13] = (int)nbparstat[0];
      for (i=0;i<nbparstat[0];i++) params[(ncolparams[0])*(row-1)+ 14 + i] = parstat[i];




    
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
    delete[] funclisttmp; 
    delete[] thetavectmp;
    delete[] xvectmp;
    delete[] ptmp;
    delete[] nptmp;
    delete[] centertmp;
    delete[] scaletmp;

    PutRNGstate();
    return Rcpp::List::create(Rcpp::Named("decision") = decision,
                          Rcpp::Named("params") = params);

}

RcppExport SEXP powcompeasyRcpp(SEXP MSEXP, SEXP paramsSEXP, SEXP ncolparamsSEXP, SEXP decisionSEXP, SEXP decisionlenSEXP,  
				SEXP modelnumSEXP, SEXP funclistSEXP, SEXP thetavecSEXP, SEXP xvecSEXP, SEXP pSEXP, SEXP npSEXP, 
				SEXP RlawsSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
  BEGIN_RCPP
  IntegerVector M = Rcpp::as<IntegerVector >(MSEXP);
  NumericVector params = Rcpp::as<NumericVector >(paramsSEXP);
  IntegerVector ncolparams = Rcpp::as<IntegerVector >(ncolparamsSEXP);
  IntegerVector decision = Rcpp::as<IntegerVector >(decisionSEXP);
  IntegerVector decisionlen = Rcpp::as<IntegerVector >(decisionlenSEXP);
  IntegerVector modelnum = Rcpp::as<IntegerVector >(modelnumSEXP);
  List funclist = Rcpp::as<List >(funclistSEXP);
  NumericVector thetavec = Rcpp::as<NumericVector >(thetavecSEXP);
  NumericVector xvec = Rcpp::as<NumericVector >(xvecSEXP);
  IntegerVector p = Rcpp::as<IntegerVector >(pSEXP);
  IntegerVector np = Rcpp::as<IntegerVector >(npSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  List Rlaws = Rcpp::as<List >(RlawsSEXP);
  List Rstats = Rcpp::as<List >(RstatsSEXP);
  SEXP __result = powcompeasyRcpp2( M,  params,  ncolparams,  decision,  decisionlen, modelnum,  funclist,  thetavec,  xvec,  p,  np,  Rlaws, Rstats, center, scale);
  return Rcpp::wrap(__result);
  END_RCPP
}




  // computation of all the statistic values in order to obtain the critical values
  SEXP compquantRcpp2(IntegerVector n, IntegerVector law, IntegerVector stat, IntegerVector M, NumericVector statvec,  
		      IntegerVector nbparlaw, NumericVector parlaw, IntegerVector nbparstat, NumericVector parstat, IntegerVector modelnum, List funclist, 
		      NumericVector thetavec, NumericVector xvec, IntegerVector p, IntegerVector np, Function Rlaw, Function Rstat, IntegerVector center, IntegerVector scale) {


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
    level[0] = 1.0;
    statistic[0] = 0.0;
	
	// We set pvalcomp[0] = 0 so the function compquant() doesn't compute p-value each time it is called
	// if you want to retrieve these p-values, please change pvalue[0] to 0
    pvalue[0] = 0.0;
    pvalcomp[0] = 0;
	
	// end
	
    critvalL[0] = 0.0;
    critvalR[0] = 0.0;

    int *usecrit, *alter;
    usecrit = new int[1];
    alter = new int[1];
    usecrit[0] = 0;
    alter[0] = 0;

    int *decisiontmp, *nblevel;
    decisiontmp = new int[1];
    decisiontmp[0] = 0;
    nblevel = new int[1];
    nblevel[0] = 1;


	     int i, j, k, *getname, *ntmp, *centertmp, *scaletmp, *nbparstattmp, *nbparlawtmp;
    getname = new int[1];
    ntmp = new int[1];
    ntmp[0] = n[0];
    centertmp = new int[1];
    centertmp[0] = center[0];
    scaletmp = new int[1];
    scaletmp[0] = scale[0];
    nbparstattmp = new int[1];
    nbparstattmp[0] = nbparstat[0];
    nbparlawtmp = new int[1];
    nbparlawtmp[0] = nbparlaw[0];

    double *x, *parstattmp, *parlawtmp;
    parstattmp = new double[nbparstat[0]];
    for (i=1;i<=nbparstat[0];i++) parstattmp[i-1] = parstattmp[i-1];
    parlawtmp = new double[nbparlaw[0]];
    for (i=1;i<=nbparlaw[0];i++) parlawtmp[i-1] = parlawtmp[i-1];


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
    GetRNGstate(); Rcpp::RNGScope __rngScope;



    getname[0] = 0;	
    if (law[0] == 0) {NumericVector mysample;
      for (i=1;i<=M[0];i++) {


        //  Rprintf("3\n ");
	  //	  Rf_PrintValue(nbparlaw[0]);

	List resultsample = gensampleRcpp2(Rlaw, n, parlaw, nbparlaw[0], "", center, scale);
	//          Rprintf("4\n ");
	//		Rf_PrintValue(resultsample["sample"]);

	   mysample = Rcpp::as<NumericVector >(resultsample["sample"]);
	  //      Rprintf("5\n ");
	//	model(modelnum[0],funclist,thetavec,xvec,n,x,p,np);   // on applique le modèle
	if (stat[0] == 0) {    
	  
	      NumericVector levelSEXP = level[0];

	      // printf("%f",level[0]);

	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstat,mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);

	  NumericVector mystatistic = resultstat["statistic"];
	  statvec[i-1] = mystatistic[0];      
	} else {  
	  x = new double[n[0]];
	  for (j=1;j<=n[0];j++) x[j-1] = mysample[j-1];
	  statcompute(stat[0], x, ntmp, level, nblevel, name2, getname, statistic, pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decisiontmp,parstattmp,nbparstattmp);
	  statvec[i-1] = statistic[0]; 
      //On libere de la memoire
	  delete[] x;
     
	}
      
      
      }
    } else {


      for (i=1;i<=M[0];i++) {
      
	x = new double[n[0]];
	for (j=1;j<=n[0];j++) x[j-1] = 0.0;
 

	gensample(law[0],ntmp,x,name1,getname,parlawtmp,nbparlawtmp,setseed,centertmp,scaletmp);
	//	model(modelnum[0],funclist,thetavec,xvec,n,x,p,np);   // on applique le modèle
     
	if (stat[0] == 0) {
	      NumericVector mysample (ntmp[0]);
	      for (k=1;k<=n[0];k++)  mysample[k-1] = x[k-1];
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstat,mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      NumericVector mystatistic = resultstat["statistic"];
	      statvec[i-1] = mystatistic[0];      
	} else { 

	  statcompute(stat[0], x, ntmp, level, nblevel, name2, getname, statistic, pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decisiontmp,parstattmp,nbparstattmp);
	  statvec[i-1] = statistic[0];      
	}
      
      //On libere de la memoire
	delete[] x;
      
      }    PutRNGstate();
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
    delete[] ntmp;
    delete[] centertmp;
    delete[] scaletmp;
    delete[] nbparstattmp;
    delete[] parstattmp;
    delete[] nbparlawtmp;
    delete[] parlawtmp;

    PutRNGstate();
    return Rcpp::List::create(Rcpp::Named("law.pars") = parlaw,
			      Rcpp::Named("nbparlaw") = nbparlaw,
			      Rcpp::Named("statvec") = statvec,
			      Rcpp::Named("stat.pars") = parstat,
			      Rcpp::Named("nbparstat") = nbparstat);
    
  }




RcppExport SEXP compquantRcpp(SEXP nSEXP, SEXP lawSEXP, SEXP statSEXP, SEXP MSEXP, SEXP statvecSEXP,  
			      SEXP nbparlawSEXP, SEXP parlawSEXP, SEXP nbparstatSEXP, SEXP parstatSEXP, SEXP modelnumSEXP, SEXP funclistSEXP, SEXP thetavecSEXP, SEXP xvecSEXP, 
			      SEXP pSEXP, SEXP npSEXP, SEXP RlawSEXP, SEXP RstatSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
  BEGIN_RCPP
  IntegerVector n = Rcpp::as<IntegerVector >(nSEXP);
  IntegerVector law = Rcpp::as<IntegerVector >(lawSEXP);
  IntegerVector stat = Rcpp::as<IntegerVector >(statSEXP);
  IntegerVector M = Rcpp::as<IntegerVector >(MSEXP);
  NumericVector statvec = Rcpp::as<NumericVector >(statvecSEXP);
  IntegerVector nbparlaw = Rcpp::as<IntegerVector >(nbparlawSEXP);
  NumericVector parlaw = Rcpp::as<NumericVector >(parlawSEXP);
  IntegerVector nbparstat = Rcpp::as<IntegerVector >(nbparstatSEXP);
  NumericVector parstat = Rcpp::as<NumericVector >(parstatSEXP);
  IntegerVector modelnum = Rcpp::as<IntegerVector >(modelnumSEXP);
  List funclist = Rcpp::as<List >(funclistSEXP);
  NumericVector thetavec = Rcpp::as<NumericVector >(thetavecSEXP);
  NumericVector xvec = Rcpp::as<NumericVector >(xvecSEXP);
  IntegerVector p = Rcpp::as<IntegerVector >(pSEXP);
  IntegerVector np = Rcpp::as<IntegerVector >(npSEXP);
  Function Rlaw = Rcpp::as<Function >(RlawSEXP);
  Function Rstat = Rcpp::as<Function >(RstatSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  SEXP __result = compquantRcpp2(n, law, stat, M, statvec, 
				 nbparlaw, parlaw, nbparstat, parstat , modelnum,  funclist,  thetavec,  xvec,  p,  np,  Rlaw, Rstat, center, scale);
  return Rcpp::wrap(__result);
  END_RCPP
}

  // Computation of the power of the test statistic
  SEXP powcompfastRcpp2(IntegerVector M, IntegerVector vectlaws, IntegerVector lawslen, IntegerVector vectn, IntegerVector vectnlen, 
			IntegerVector vectstats, IntegerVector statslen,IntegerVector  decision,IntegerVector  decisionlen, NumericVector level,
			IntegerVector nblevel,NumericVector  critvalL,NumericVector  critvalR, IntegerVector usecrit, IntegerVector alter,
			IntegerVector nbparlaw, NumericVector parlaw,IntegerVector  nbparstat,NumericVector  parstat, IntegerVector modelnum, List funclist, 
			NumericVector thetavec, NumericVector xvec,IntegerVector  p,IntegerVector  np, List Rlaws, List Rstats, IntegerVector center, 
			IntegerVector scale, IntegerVector compquant) {
	      
    // Warning: When compquant[0] == 1, critvalL should be initialized with lawslen[0] * M[0] * vectnlen[0] * statslen[0] double values since it will contain (when output)
    // all the test statistics generated. Most probably, lawslen[0] should be equal to 1.

    int gensample(int law, int *xlen, double *x, char **name1, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name2, int *getname, double *statistic, 
		     int *pvalcomp,double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);
    int model(int modelnum, char** funclist, double *thetavec, double *xvec, int *xlen, double *x, int *p, int *np);
   
    double *statistic, *pvalue; // POUR L'INSTANT JE N'EN FAIT RIEN de pvalue!! Si je veux les récupérer dans R il faudra faire des modifs!! A voir ...
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
    }    
    
    int *decisiontmp;
    decisiontmp = new int[nblevel[0]];
    for (i = 0; i < nblevel[0]; i++) decisiontmp[i] = 0;
    
    int *setseed;
    setseed = new int[1];
    setseed[0] = 0;
    GetRNGstate(); Rcpp::RNGScope __rngScope;
	      
	      
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
    int kmax=0; 		// kmax = parstats.len.max in powcomp-fast.R
    for (t = 0; t < statslen[0]; t++) {
      if (kmax <= nbparstat[t]) {
	kmax = nbparstat[t];
      } //else kmax = kmax;
    }
    if (kmax == 0) kmax = 1;
    parstattmp = new double[kmax];
	      
    int stlen1, stlen2;
	      
    getname[0] = 0;	

    double *leveltmp;
    int ii, k, *nbleveltmp, *centertmp, *scaletmp;
    nbleveltmp = new int[1];
    nbleveltmp[0] = nblevel[0];
    centertmp = new int[1];
    centertmp[0] = center[0];
    scaletmp = new int[1];
    scaletmp[0] = scale[0];
    leveltmp = new double[nblevel[0]];
    for (i = 0; i < nblevel[0]; i++) leveltmp[i] = level[i];

    for (law = 0; law < lawslen[0]; law++) {

      for (i = 1; i <= M[0]; i++) {	      
	      
	xlen[0] = maxn; // on génère un échantillon de taille maxn de loi law 
	      
	nbparlawtmp[0] = nbparlaw[law];
	parlawtmp[0] = parlaw[0 + 4 * law];
	parlawtmp[1] = parlaw[1 + 4 * law];
	parlawtmp[2] = parlaw[2 + 4 * law];
	parlawtmp[3] = parlaw[3 + 4 * law];
	      
	NumericVector mysample(xlen[0]);
	if (vectlaws[law] == 0) {
	  Function rlawfunc = Rcpp::as<Function >(Rlaws[law]);
	  IntegerVector nn(1);
	  nn[0] = maxn;
	  NumericVector parlawbis(4);
	  parlawbis[0] = parlaw[0 + 4 * law];
	  parlawbis[1] = parlaw[1 + 4 * law];
	  parlawbis[2] = parlaw[2 + 4 * law];
	  parlawbis[3] = parlaw[3 + 4 * law];

	  List resultsample = gensampleRcpp2(rlawfunc, nn, parlawbis, nbparlaw[law], "",center,scale);
	  mysample = resultsample["sample"];
	  for (j = 0; j < maxn; j++) x[j] = mysample[j];
	      //	model(modelnum[0],funclist,thetavec,xvec,xlen,x,p,np);   // on applique le modèle
	} else {
	  gensample(vectlaws[law], xlen, x, name1, getname, parlawtmp, nbparlawtmp, setseed, centertmp, scaletmp);
	  for (k = 0; k < xlen[0]; k++) mysample[k] = x[k];
	      
	      //	  model(modelnum[0],funclist,thetavec,xvec,xlen,x,p,np);   // on applique le modèle
	  PutRNGstate();
	}
	      
	if (i==1) {
	  nbparlaw[law] = nbparlawtmp[0];
	  parlaw[0 + 4 * law] = parlawtmp[0];
	  parlaw[1 + 4 * law] = parlawtmp[1];
	  parlaw[2 + 4 * law] = parlawtmp[2];
	  parlaw[3 + 4 * law] = parlawtmp[3];
	}
	      
	
	NumericVector mydecision(nblevel[0]);
	
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
	      
	    if (vectstats[stat] == 0) {
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrittmp[0];
	      NumericVector critvalLSEXP(nblevel[0]);
	      NumericVector critvalRSEXP(nblevel[0]);
	      NumericVector mysamplebis(xlen[0]);
	      for (ii = 0; ii < xlen[0]; ii++) mysamplebis[ii] = mysample[ii];
	      for (ii = 0; ii < nblevel[0]; ii++) {
		critvalLSEXP[ii] = critvalLtmp[ii];
		critvalRSEXP[ii] = critvalRtmp[ii];
	      }
	      List resultstat = statcomputeRcpp2(Rstats[stat],mysamplebis,level,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      mydecision = resultstat["decision"];
	      if (compquant[0] == 1) {
		critvalL[i + M[0] * n + M[0] * vectnlen[0] * stat - 1] = resultstat["statistic"];
	      }
	    } else {
	      // decisiontmp est de longueur nblevel[0]
	      statcompute(vectstats[stat], x, xlen, leveltmp, nbleveltmp, name2, getname, statistic, pvalcomp, pvalue, critvalLtmp, critvalRtmp, usecrittmp, altertmp, decisiontmp, parstattmp, nbparstattmp);
	      if (compquant[0] == 1) {
		critvalL[i + M[0] * n + M[0] * vectnlen[0] * stat - 1] = statistic[0];
	      }
	      for (k = 0; k < nblevel[0]; k++) mydecision[k] = decisiontmp[k];
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
		if (vectstats[stat] == 0) {
		  decision[indtmp] = decision[indtmp] + mydecision[j];
		} else {
		  decision[indtmp] = decision[indtmp] + decisiontmp[j];
		}
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
    delete[] nbleveltmp;
    delete[] leveltmp;
    delete[] centertmp;
    delete[] scaletmp;
    
    PutRNGstate();
    return Rcpp::List::create(Rcpp::Named("M") = M,
			      Rcpp::Named("law.indices") = vectlaws,
			      Rcpp::Named("vectn") = vectn,
			      Rcpp::Named("stat.indices") = vectstats,
			      Rcpp::Named("decision") = decision,
			      Rcpp::Named("levels") = level,
			      Rcpp::Named("cL") = critvalL,
			      Rcpp::Named("cR") = critvalR,
			      Rcpp::Named("usecrit") = usecrit,
			      Rcpp::Named("alter") = alter,
			      Rcpp::Named("nbparlaws") = nbparlaw,
			      Rcpp::Named("parlaws") = parlaw,
			      Rcpp::Named("nbparstats") = nbparstat);
  }
 



RcppExport SEXP powcompfastRcpp(SEXP MSEXP, SEXP vectlawsSEXP, SEXP lawslenSEXP, SEXP vectnSEXP, SEXP vectnlenSEXP,  
				SEXP vectstatsSEXP, SEXP statslenSEXP, SEXP decisionSEXP, SEXP decisionlenSEXP, SEXP levelSEXP, SEXP nblevelSEXP, SEXP critvalLSEXP, SEXP critvalRSEXP, 
				SEXP usecritSEXP, SEXP alterSEXP, SEXP nbparlawSEXP, SEXP parlawSEXP, SEXP nbparstatSEXP, SEXP parstatSEXP, SEXP modelnumSEXP, SEXP funclistSEXP, 
				SEXP thetavecSEXP, SEXP xvecSEXP, SEXP pSEXP, SEXP npSEXP, SEXP RlawsSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP compquantSEXP) {
  BEGIN_RCPP
  IntegerVector M = Rcpp::as<IntegerVector >(MSEXP);
  IntegerVector vectlaws = Rcpp::as<IntegerVector >(vectlawsSEXP);
  IntegerVector lawslen = Rcpp::as<IntegerVector >(lawslenSEXP);
  IntegerVector vectn = Rcpp::as<IntegerVector >(vectnSEXP);
  IntegerVector vectnlen = Rcpp::as<IntegerVector >(vectnlenSEXP);
  IntegerVector vectstats = Rcpp::as<IntegerVector >(vectstatsSEXP);
  IntegerVector statslen = Rcpp::as<IntegerVector >(statslenSEXP);
  IntegerVector decision = Rcpp::as<IntegerVector >(decisionSEXP);
  IntegerVector decisionlen = Rcpp::as<IntegerVector >(decisionlenSEXP);
  NumericVector level = Rcpp::as<NumericVector >(levelSEXP);
  IntegerVector nblevel = Rcpp::as<IntegerVector >(nblevelSEXP);
  NumericVector critvalL = Rcpp::as<NumericVector >(critvalLSEXP);
  NumericVector critvalR = Rcpp::as<NumericVector >(critvalRSEXP);
  IntegerVector usecrit = Rcpp::as<IntegerVector >(usecritSEXP);
  IntegerVector alter = Rcpp::as<IntegerVector >(alterSEXP);
  IntegerVector nbparlaw = Rcpp::as<IntegerVector >(nbparlawSEXP);
  NumericVector parlaw = Rcpp::as<NumericVector >(parlawSEXP);
  IntegerVector nbparstat = Rcpp::as<IntegerVector >(nbparstatSEXP);
  NumericVector parstat = Rcpp::as<NumericVector >(parstatSEXP);
  IntegerVector modelnum = Rcpp::as<IntegerVector >(modelnumSEXP);
  List funclist = Rcpp::as<List >(funclistSEXP);
  NumericVector thetavec = Rcpp::as<NumericVector >(thetavecSEXP);
  NumericVector xvec = Rcpp::as<NumericVector >(xvecSEXP);
  IntegerVector p = Rcpp::as<IntegerVector >(pSEXP);
  IntegerVector np = Rcpp::as<IntegerVector >(npSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  IntegerVector compquant = Rcpp::as<IntegerVector >(compquantSEXP);
  List Rlaws = Rcpp::as<List >(RlawsSEXP);
  List Rstats = Rcpp::as<List >(RstatsSEXP);
  SEXP __result = powcompfastRcpp2(M, vectlaws, lawslen, vectn, vectnlen, vectstats, statslen, decision, decisionlen, level, nblevel,
				   critvalL, critvalR, usecrit, alter, nbparlaw, parlaw, nbparstat, parstat, modelnum, funclist, thetavec, 
				   xvec, p, np, Rlaws, Rstats, center, scale, compquant);
  return Rcpp::wrap(__result);
  END_RCPP
}


  // Computation of the p-values matrix
  SEXP matrixpvalRcpp2(IntegerVector N,  IntegerVector lawindex,  IntegerVector xlen,  IntegerVector nbparams,  NumericVector parlaw,  IntegerVector statindices,  
		       IntegerVector nbstats, IntegerVector altervec, NumericVector parstatmultvec, IntegerVector nbparstatvec, NumericVector res, Function Rlaw, List Rstats, IntegerVector center, IntegerVector scale) {

    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name2, int *getname, double *statistic, 
		     int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);

    int i, j, k;

    int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);

    double *x; 
    x = new double[xlen[0]];
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

    int *getname, *centertmp, *scaletmp;
    getname = new int[1];
    getname[0] = 0;
    centertmp = new int[1];
    centertmp[0] = center[0];
    scaletmp = new int[1];
    scaletmp[0] = scale[0];

    int *setseed;
    setseed = new int[1];
    setseed[0] = 1;
    GetRNGstate(); Rcpp::RNGScope __rngScope;

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

    int * nbparamstat, *xlentmp, *nbparamstmp;
    nbparamstat = new int[1];
    xlentmp = new int[1];
    xlentmp[0] = xlen[0];
    nbparamstmp = new int[1];
    nbparamstmp[0] = nbparams[0];

    IntegerVector nSEXP = xlen[0]; 
    NumericVector mysample(xlen[0]);

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

      if (lawindex[0] == 0) {Rcpp::RNGScope __rngScope;
	List resultsample = gensampleRcpp2(Rlaw, nSEXP, parlaw, nbparams[0], "",center,scale);
	mysample = resultsample["sample"];
	for (k=1;k<=xlen[0];k++) x[k-1] = mysample[k-1];

      } else {
	GetRNGstate(); 
	gensample(lawindex[0],xlentmp,x,name,getname,params,nbparamstmp,setseed,centertmp,scaletmp);  // We retrieve x
	for (k=1;k<=xlen[0];k++) mysample[k-1] = x[k-1];
	PutRNGstate();
      }
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
	  if (statindex == 0) {
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[i],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      pvalue[0] = resultstat["pvalue"];
	    } else {
	      statcompute(statindex, x,xlentmp,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparamstat);
	    }
	  delete[] paramstat;
	} else {

	      if (statindex == 0) {
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[i],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      pvalue[0] = resultstat["pvalue"];
	    } else {
	      statcompute(statindex,x,xlentmp,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,(double*)0,nbparamstat);
	    }
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
    delete[] xlentmp;
    delete[] centertmp;
    delete[] scaletmp;

    PutRNGstate();

    return res;

  }
 

RcppExport SEXP matrixpvalRcpp(SEXP NSEXP, SEXP lawindexSEXP, SEXP xlenSEXP, SEXP nbparamsSEXP, SEXP parlawSEXP,  
			       SEXP statindicesSEXP, SEXP nbstatsSEXP, SEXP altervecSEXP, SEXP parstatmultvecSEXP, SEXP nbparstatvecSEXP, SEXP resSEXP, SEXP RlawSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
  BEGIN_RCPP
  IntegerVector N = Rcpp::as<IntegerVector >(NSEXP);
  IntegerVector lawindex = Rcpp::as<IntegerVector >(lawindexSEXP);
  IntegerVector xlen = Rcpp::as<IntegerVector >(xlenSEXP);
  IntegerVector nbparams = Rcpp::as<IntegerVector >(nbparamsSEXP);
  NumericVector parlaw = Rcpp::as<NumericVector >(parlawSEXP);
  IntegerVector statindices = Rcpp::as<IntegerVector >(statindicesSEXP);
  IntegerVector nbstats = Rcpp::as<IntegerVector >(nbstatsSEXP);
  IntegerVector altervec = Rcpp::as<IntegerVector >(altervecSEXP);
  NumericVector parstatmultvec = Rcpp::as<NumericVector >(parstatmultvecSEXP);
  IntegerVector nbparstatvec = Rcpp::as<IntegerVector >(nbparstatvecSEXP);
  NumericVector res = Rcpp::as<NumericVector >(resSEXP);
  Function Rlaw = Rcpp::as<Function >(RlawSEXP);
  List Rstats = Rcpp::as<List >(RstatsSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  SEXP __result = matrixpvalRcpp2( N,  lawindex,  xlen,  nbparams,  parlaw,  statindices,  nbstats,  altervec,  parstatmultvec,  nbparstatvec,  res, Rlaw, Rstats, center, scale);
  return Rcpp::wrap(__result);
  END_RCPP
}


  // Computation of the p-values matrix using Monte-Carlo
  SEXP matrixpvalMCRcpp2(IntegerVector n, IntegerVector lawindex, IntegerVector nbstats,IntegerVector  M, IntegerVector statindices,  
				IntegerVector  nbparstatvec, NumericVector parstatmultvec, List funclist, IntegerVector N, 
				IntegerVector  nulldist, IntegerVector nbparams, IntegerVector altervec, NumericVector parstat, 
			 IntegerVector  nbparstat, NumericVector res, Function Rlawindex, Function Rnulldist, List Rstats, IntegerVector center, IntegerVector scale) {


    void statcompute(int stat, double *x, int *xlen, double *level, int *nblevel, char **name2, int *getname, double *statistic, 
		     int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat);

    void compquantc(int *n, int *law, int *stat, int *M, double *statvec,// char **lawname, char **statname, 
    		   int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale);
    int gensample(int law, int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed, int *center, int *scale);
    int i, j, k;

    int *stat, *ntmp, *lawindextmp, *Mtmp, *nbparamstmp, *centertmp, *scaletmp;
    stat = new int[1];
    stat[0] = 0;
    ntmp = new int[1];
    ntmp[0] = n[0];
    lawindextmp = new int[1];
    lawindextmp[0] = lawindex[0];
    Mtmp = new int[1];
    Mtmp[0] = M[0];
    nbparamstmp = new int[1];
    nbparamstmp[0] = nbparams[0];
    centertmp = new int[1];
    centertmp[0] = center[0];
    scaletmp = new int[1];
    scaletmp[0] = scale[0];

    double *statvec;
    statvec = new double[M[0]];
    for (i=0;i<M[0];i++) statvec[i] = 0.0;


    char **lawname;         // lawname=paste(rep(" ",50), sep = "", collapse = "")
    lawname = new char*[1];
    lawname[0] = new char[50];
    for (i=0;i<50;i++) lawname[0][i] = ' ';

    char **statname;          // statname = paste(rep(" ",50), sep = "", collapse = "")
    statname = new char*[1];
    statname[0] = new char[50];
    for (i=0;i<50;i++) statname[0][i] = ' ';

    char **funclisttmp;
    funclisttmp = new char*[1];
    funclisttmp[0] = (char *)R_ExternalPtrAddr(funclist[0]);

    int *nbparlaw;
    nbparlaw = new int[1];
    nbparlaw[0] = 0;

    double *parlaw;
    parlaw = new double[4];
    parlaw[0] = 0.0; parlaw[1] = 0.0; parlaw[2] = 0.0; parlaw[3] = 0.0; 

    NumericVector parlawSEXP(4);
    parlawSEXP[0] = 0.0; parlawSEXP[1] = 0.0; parlawSEXP[2] = 0.0; parlawSEXP[3] = 0.0; 

    int *nbparamstat;
    nbparamstat = new int[1];
    nbparamstat[0] = 0;

    int cmpt;

    int *modelnum;
    modelnum = new int[1];
    modelnum[0] = 1;

    IntegerVector modelnumSEXP(1);
    modelnumSEXP[0] = 1;

    double *thetavec;
    thetavec = new double[1];
    thetavec[0] = 0.0;

    NumericVector thetavecSEXP = 0.0;

    double *xvec;
    xvec = new double[1];
    xvec[0] = 0.0;

    NumericVector xvecSEXP = 0.0;

    int *p;
    p = new int[1];
    p[0] = 1;

    IntegerVector pSEXP(1);
    pSEXP[0] = 1;

    int *np;
    np = new int[1];
    np[0] = 1;

    IntegerVector npSEXP(1);
    npSEXP = 1;

    double *liststat;
    liststat = new double[M[0]*nbstats[0]];

      cmpt = 0;

      IntegerVector statSEXP(1);
    NumericVector statvecSEXP(M[0]);
 
    for (i=0;i<nbstats[0];i++) {

      stat[0] = statindices[i];
      Function Rstat = Rstats[i];
      statSEXP[0] = stat[0];
   
      for (j=0;j<M[0];j++) statvecSEXP[j] = 0.0;

      nbparamstat[0] = nbparstatvec[i];

      if (nbparamstat[0]>0) {	    
	double * paramstat;
	paramstat = new double[nbparamstat[0]];    
	for (k=0;k<nbparamstat[0];k++) paramstat[k] = parstatmultvec[cmpt+k];

	cmpt = cmpt + nbparamstat[0];
	if (lawindex[0] == 0) {
	  compquantRcpp2(n, lawindex, statSEXP, M, statvecSEXP, 
			 nbparams, parlawSEXP, nbparstat, parstat , modelnumSEXP,  funclist,  thetavecSEXP,  xvecSEXP,  pSEXP,  npSEXP,  Rlawindex, Rstat, center, scale);
	} else {    
	  compquantc(ntmp,lawindextmp,stat,Mtmp,statvec,//lawname,statname,
		    nbparlaw,parlaw,nbparamstat,paramstat,modelnum,funclisttmp,thetavec,xvec,p,np, centertmp, scaletmp);
	}
	delete[] paramstat;
      } else {    
	if (lawindex[0] == 0) {

	  compquantRcpp2(n, lawindex, statSEXP, M, statvecSEXP, 
			 nbparams, parlawSEXP, nbparstat, parstat , modelnumSEXP,  funclist,  thetavecSEXP,  xvecSEXP,  pSEXP,  npSEXP,  Rlawindex, Rstat, center, scale);
	} else {
	  
	  compquantc(ntmp,lawindextmp,stat,Mtmp,statvec,// lawname,statname,
		    nbparlaw,parlaw,nbparamstat,(double*)0,modelnum,funclisttmp,thetavec,xvec,p,np, centertmp, scaletmp);
	}
      }

      for (j=0;j<M[0];j++) liststat[i*M[0] + j] = statvec[j];
   
    }




    double *x; 
    x = new double[n[0]];
    for (i=0;i<n[0];i++) x[i] = 0.0; 


    double *params; // ATTENTION. Bizarre!! NE semble pas vraiment utile car parlaw a ete initialise plus haut a 0 pour toutes ses composantes!!!
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

    NumericVector paramsSEXP(4);
    if (nbparams[0] == 0) {
      paramsSEXP[0] = 0.0;
      paramsSEXP[1] = 0.0;
      paramsSEXP[2] = 0.0;
      paramsSEXP[3] = 0.0;
    }
    if (nbparams[0] == 1) {
      paramsSEXP[0] = parlaw[0];
      paramsSEXP[1] = 0.0;
      paramsSEXP[2] = 0.0;
      paramsSEXP[3] = 0.0;
    }
    if (nbparams[0] == 2) {
      paramsSEXP[0] = parlaw[0];
      paramsSEXP[1] = parlaw[1];
      paramsSEXP[2] = 0.0;
      paramsSEXP[3] = 0.0;
    }
    if (nbparams[0] == 3) {
      paramsSEXP[0] = parlaw[0];
      paramsSEXP[1] = parlaw[1];
      paramsSEXP[2] = parlaw[2];
      paramsSEXP[3] = 0.0;
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
    GetRNGstate(); Rcpp::RNGScope __rngScope;


    int *nblevel, *usecrit, *decision;
    nblevel = new int[1];
    nblevel[0] = 0;
    usecrit = new int[1];
    usecrit[0] = 0;
    decision = new int[1];
    decision[0] = 0;
    double *level, *critvalL, *critvalR, *statistic, *pvalue;
    int *pvalcomp;
    level = new double[1];
    level[0] = 0.0;
    critvalL = new double[1];
    critvalL[0] = 0.0;
    critvalR = new double[1];
    critvalR[0] = 0.0;
    statistic = new double[1];
    statistic[0] = 0.0;
    pvalue = new double[1];
    pvalue[0] = 0.0;
    pvalcomp = new int[1];
    pvalcomp[0] = 0;


    double *liststati;
    liststati = new double[M[0]];
    for (k=0;k<M[0];k++) liststati[k] = 0.0;

    double q2;

    double meanstatsup, meanstatinf;

    int *alter;
    alter = new int[1];
    alter[0] = 0;

    NumericVector mysample(n[0]);


    for (j=0;j<N[0];j++) {

      if (nulldist[0] == 0) {Rcpp::RNGScope __rngScope;
	List resultsample = gensampleRcpp2(Rnulldist, n, paramsSEXP, nbparlaw[0], "",center,scale);
	mysample = resultsample["sample"];
	for (k=1;k<=n[0];k++) x[k-1] = mysample[k-1];
      } else {
	GetRNGstate(); 
	gensample(nulldist[0],ntmp,x,name,getname,params,nbparamstmp,setseed,centertmp,scaletmp);  // We retrieve x
	for (k=1;k<=n[0];k++)  mysample[k-1] = x[k-1];
	PutRNGstate();
      }

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
	  if (statindices[i] == 0) {
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[i],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      statistic[0] = resultstat["statistic"];
	      pvalue[0] = resultstat["pvalue"];
	    } else {
	      statcompute(statindices[i],x,ntmp,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,paramstat,nbparamstat);    // We obtain one test statistic
	    }
	  delete[] paramstat;
	} else {
	  if (statindices[i] == 0) {
	      NumericVector levelSEXP = level[0];
	      IntegerVector usecritSEXP(1);
	      usecritSEXP[0] = usecrit[0];
	      NumericVector critvalLSEXP = critvalL[0];
	      NumericVector critvalRSEXP = critvalR[0];
	      List resultstat = statcomputeRcpp2(Rstats[i],mysample,levelSEXP,usecritSEXP,critvalLSEXP,critvalRSEXP);
	      statistic[0] = resultstat["statistic"];
	      pvalue[0] = resultstat["pvalue"];
	    } else {
	      statcompute(statindices[i],x,ntmp,level,nblevel,name,getname,statistic,pvalcomp,pvalue,critvalL,critvalR,usecrit,alter,decision,(double*)0,nbparamstat);    // We obtain one test statistic
	    }
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
      delete[] name[i-1];
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
    delete[] ntmp;
    delete[] lawindextmp;
    delete[] Mtmp;
    delete[] centertmp;
    delete[] scaletmp;
    delete[] nbparamstmp;
    delete[] funclisttmp;

    PutRNGstate();
    return res;

  }



RcppExport SEXP matrixpvalMCRcpp(SEXP nSEXP, SEXP lawindexSEXP, SEXP nbstatsSEXP, SEXP MSEXP, SEXP statindicesSEXP,  
				 SEXP nbparstatvecSEXP, SEXP parstatmultvecSEXP, SEXP funclistSEXP, SEXP NSEXP, 
				 SEXP nulldistSEXP, SEXP nbparamsSEXP, SEXP altervecSEXP, SEXP parstatSEXP, 
				 SEXP nbparstatSEXP, SEXP resSEXP, SEXP RlawindexSEXP, SEXP RnulldistSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP) {
  BEGIN_RCPP
  IntegerVector n = Rcpp::as<IntegerVector >(nSEXP);
  IntegerVector lawindex = Rcpp::as<IntegerVector >(lawindexSEXP);
  IntegerVector nbstats = Rcpp::as<IntegerVector >(nbstatsSEXP);
  IntegerVector M = Rcpp::as<IntegerVector >(MSEXP);
  IntegerVector statindices = Rcpp::as<IntegerVector >(statindicesSEXP);
  IntegerVector nbparstatvec = Rcpp::as<IntegerVector >(nbparstatvecSEXP);
  NumericVector parstatmultvec = Rcpp::as<NumericVector >(parstatmultvecSEXP);
  List funclist = Rcpp::as<List >(funclistSEXP);
  IntegerVector N = Rcpp::as<IntegerVector >(NSEXP);
  IntegerVector nulldist = Rcpp::as<IntegerVector >(nulldistSEXP);
  IntegerVector nbparams = Rcpp::as<IntegerVector >(nbparamsSEXP);
  IntegerVector altervec = Rcpp::as<IntegerVector >(altervecSEXP);
  NumericVector parstat = Rcpp::as<NumericVector >(parstatSEXP);
  IntegerVector nbparstat = Rcpp::as<IntegerVector >(nbparstatSEXP);
  NumericVector res = Rcpp::as<NumericVector >(resSEXP);
  Function Rlawindex = Rcpp::as<Function >(RlawindexSEXP);
  Function Rnulldist = Rcpp::as<Function >(RnulldistSEXP);
  List Rstats = Rcpp::as<List >(RstatsSEXP);
  IntegerVector center = Rcpp::as<IntegerVector >(centerSEXP);
  IntegerVector scale = Rcpp::as<IntegerVector >(scaleSEXP);
  SEXP __result = matrixpvalMCRcpp2( n,  lawindex,  nbstats,  M,  statindices,  
				  nbparstatvec,  parstatmultvec,  funclist,  N, 
				  nulldist,  nbparams,  altervec,  parstat, 
				     nbparstat,  res,  Rlawindex, Rnulldist, Rstats, center, scale );
  return Rcpp::wrap(__result);
  END_RCPP
}




  
}
