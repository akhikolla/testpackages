/*
cd PoweR/
R-devel
> tools::package_native_routine_registration_skeleton(".")
*/
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/*
  The following symbols/expressions for .NAME have been omitted

    dontCheck(stat.name)
    dontCheck(Claw.name)
    dontCheck(Cstat.name)
    paste("stat", stat.index, sep = "")

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void calcfx(double *pvals, int *pvalslen, double *xi, int *xilen, double *fx);
extern void compquantc(int *n, int *law, int *stat, int *M, double *statvec, //char **lawname, char **statname, 
		 int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale);
extern void matrixpval(int* N, int* lawindex, int *xlen, int *nbparams, double *parlaw, int *statindices, int* nbstats, int *altervec, double *parstatmultvec, int *nbparstatvec, double *res, int *center, int *scale);
extern void matrixpvalMC(int *n, int *lawindex, int* nbstats, int *M, int *statindices, int *nbparstatvec, double *parstatmultvec, char** funclist, int *N, int *nulldist, int *nbparams, int *altervec, double *parstat, int *nbparstat, double *res, int *center, int *scale);
extern void powcompeasy(int *M, double *params, int *ncolparams, int *decision, int *decisionlen, //char **lawnames, char **statnames, 
		   int *modelnum, char** funclist, double *thetavec, double *xvec, int *p, int *np, int *center, int *scale);
extern void powcompfast(int *M, int *vectlaws, int *lawslen, int *vectn, int *vectnlen, int *vectstats, int *statslen, int *decision, int *decisionlen, double *level, int *nblevel,
		   double *critvalL, double *critvalR, int *usecrit, int *alter, int *nbparlaw, double *parlaw, int *nbparstat, double *parstat, int *modelnum, char** funclist, 
		   double *thetavec, double *xvec, int *p, int *np, int *center, int *scale, int *compquant);

/* .Call calls */
extern SEXP compquantRcpp(SEXP nSEXP, SEXP lawSEXP, SEXP statSEXP, SEXP MSEXP, SEXP statvecSEXP,  
			      SEXP nbparlawSEXP, SEXP parlawSEXP, SEXP nbparstatSEXP, SEXP parstatSEXP, SEXP modelnumSEXP, SEXP funclistSEXP, SEXP thetavecSEXP, SEXP xvecSEXP, 
			      SEXP pSEXP, SEXP npSEXP, SEXP RlawSEXP, SEXP RstatSEXP, SEXP centerSEXP, SEXP scaleSEXP);
extern SEXP gensampleRcpp(SEXP rlawfuncSEXP, SEXP nSEXP, SEXP paramsSEXP , SEXP nbparamsSEXP, SEXP lawnameSEXP, SEXP centerSEXP, SEXP scaleSEXP);
extern SEXP matrixpvalMCRcpp(SEXP nSEXP, SEXP lawindexSEXP, SEXP nbstatsSEXP, SEXP MSEXP, SEXP statindicesSEXP,  
				 SEXP nbparstatvecSEXP, SEXP parstatmultvecSEXP, SEXP funclistSEXP, SEXP NSEXP, 
				 SEXP nulldistSEXP, SEXP nbparamsSEXP, SEXP altervecSEXP, SEXP parstatSEXP, 
				 SEXP nbparstatSEXP, SEXP resSEXP, SEXP RlawindexSEXP, SEXP RnulldistSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP);
extern SEXP matrixpvalRcpp(SEXP NSEXP, SEXP lawindexSEXP, SEXP xlenSEXP, SEXP nbparamsSEXP, SEXP parlawSEXP,  
			       SEXP statindicesSEXP, SEXP nbstatsSEXP, SEXP altervecSEXP, SEXP parstatmultvecSEXP, SEXP nbparstatvecSEXP, SEXP resSEXP, SEXP RlawSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP);
extern SEXP powcompeasyRcpp(SEXP MSEXP, SEXP paramsSEXP, SEXP ncolparamsSEXP, SEXP decisionSEXP, SEXP decisionlenSEXP,  
				SEXP modelnumSEXP, SEXP funclistSEXP, SEXP thetavecSEXP, SEXP xvecSEXP, SEXP pSEXP, SEXP npSEXP, 
				SEXP RlawsSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP);
extern SEXP powcompfastRcpp(SEXP MSEXP, SEXP vectlawsSEXP, SEXP lawslenSEXP, SEXP vectnSEXP, SEXP vectnlenSEXP,  
				SEXP vectstatsSEXP, SEXP statslenSEXP, SEXP decisionSEXP, SEXP decisionlenSEXP, SEXP levelSEXP, SEXP nblevelSEXP, SEXP critvalLSEXP, SEXP critvalRSEXP, 
				SEXP usecritSEXP, SEXP alterSEXP, SEXP nbparlawSEXP, SEXP parlawSEXP, SEXP nbparstatSEXP, SEXP parstatSEXP, SEXP modelnumSEXP, SEXP funclistSEXP, 
				SEXP thetavecSEXP, SEXP xvecSEXP, SEXP pSEXP, SEXP npSEXP, SEXP RlawsSEXP, SEXP RstatsSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP compquantSEXP);
extern SEXP statcomputeRcpp(SEXP rstatfuncSEXP, SEXP echSEXP, SEXP levelsSEXP, SEXP usecritSEXP, SEXP critvalLSEXP, SEXP critvalRSEXP);

