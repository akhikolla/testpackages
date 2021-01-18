#include "Krigtypes.h"
#define R_NO_REMAP
#include "smooth.h" // (includes qr.h which includes R.h)
#include "smoothFriends.h"
#include "R.h" // fn error() pour message erreur
using namespace Rcpp;


std::vector< CSmooth* >CKrigptrTable(0);
int fittedparamnbr=0; // will be better filled later
namespace NS_GG {
int a=3;
covTypedef b; // found no way to make templates here
}

// [[Rcpp::export]]
SEXP newCSmooth( // ici la def
    SEXP xy,
    SEXP nrowxy, //with replicates
    SEXP ncolxy,
    SEXP nuniquerows, // required to tell the allocated size of c, d, D, u in *R*
    SEXP GCV,
    SEXP optimiseBool,
    SEXP verbosity
) {
  NumericVector xxy(xy);
  int xnrowxy=as<int>(nrowxy);
  int xncolxy=as<int>(ncolxy);
  int xnuniquerows=as<int>(nuniquerows);
  double xGCV=as<double>(GCV);
  int xoptimiseBool=as<int>(optimiseBool);
  int xverbosity=as<int>(verbosity);
  bool success=intern_newCSmooth( //
    xxy.begin(),
    &xnrowxy, //with replicates
    &xncolxy,
    &xnuniquerows,
    &xGCV,
    &xoptimiseBool,
    &xverbosity
  );
  return(wrap(success));
}

bool intern_newCSmooth( // ici la def
  double *xy,
  int *nrowxy, //with replicates
  int *ncolxy,
  int *nuniquerows, // required to tell the allocated size of c, d, D, u in *R*
  //    double *maxSmoothness,
  //    double *c,
  //    double *d,
  //    double *D,
  //    double *u,
  //    double *lambda,
  double *GCV,
  //    double *covfnparam, //p+1 form
  // int *fneval,
  int *optimiseBool,
  int *verbosity
  //    double *initcovfnparam //p+1 form  // new arg 2015/10/24
) {
  fittedparamnbr=(*ncolxy)-1;
  fnevalcounter=0;
  Cpointls local(xy,nrowxy,ncolxy);
  if((*optimiseBool) ) {
    if(*GCV==0) {
#ifdef NO_R_CONSOLE
      std::cout<<"Estimating missing parameters via generalized cross-validation...";
#else
      if (*verbosity) REprintf("Estimating missing parameters via generalized cross-validation...\n");
#endif
    } else {
#ifdef NO_R_CONSOLE
      std::cout<<"Estimating missing parameters via match of MSE estimates...";
#else
      REprintf("Estimating missing parameters via match of MSE estimates...\n");
#endif
    }
  }
  test=new CSmooth(local,*GCV,*verbosity); //constructor performs QR decomp of Tmatrix etc.
  // note that the delete is suppressed at the end of this function. cf CdeleteCKrig(...)
  int Cuniquerows=test->xy.size(); // since sort_compress() has been called by the constructor
  bool success= ((*nuniquerows)==Cuniquerows);
  if (!success) {
#ifdef NO_R_CONSOLE
    // if close impossible
#else
    std::stringstream stst;
    stst<<"(!) From intern_newCSmooth() in DLL: C code counted "<<Cuniquerows<<" unique coordinates while R declared "<<*nuniquerows<<" ones ('nuniquerows' argument)\n";
    //Rf_error(stst.str().c_str());
    REprintf(stst.str().c_str());
    REprintf("This has occurred in normal usage (as R and C algorithms for counting unique values differ)\n");
    REprintf(" but might also indicate wrong input from R (although this has never occurred).\n");
#endif
  }
  return(success);
}

// [[Rcpp::export]]
int deleteCSmooth() {// useful for deleting a pointer not stored in CKrigptrTable
  delete test;
  return(0);
}

// [[Rcpp::export]]
int flushCSmoothTable() { // useful for deleting pointers stored in CKrigptrTable
  for (std::vector< CSmooth* >::iterator it=CKrigptrTable.begin();it != CKrigptrTable.end();it++) delete (*it);
  CKrigptrTable.resize(0);
  return(0);
}



// [[Rcpp::export]]
SEXP GCV_lamVar_covFix_Wrapper( SEXP a, SEXP fixedSmoothness, SEXP returnFnvalue ){ // ici la DEF
  Rcpp::NumericVector xa(a);
  Rcpp::NumericVector xfixedSmoothness(fixedSmoothness); // may be of zero length
  bool xreturnFnvalue=as<bool>(returnFnvalue);
  double resu;
  std::vector<double> xv(0);
  for (Rcpp::NumericVector::iterator ii=xa.begin();ii!=xa.end();ii++) xv.push_back(*ii);
  if (xfixedSmoothness.size()==1) xv.push_back(xfixedSmoothness[0]); // else smoothness is already in vector
  double fnvalue=test->GCV_lamVar_covFix<double,double>(xv);
  if (xreturnFnvalue) {
    resu = fnvalue;    // GCV value at lambdaEst for fixed covparam
    if (test->verbosity && (fnevalcounter % xv.size()==0)) {
#ifdef NO_R_CONSOLE
    std::cout<<std::endl;
    for (int i=0;i<int(xv.size());i++) {std::cout<<xv[i]<<" ";}
    std::cout<<" f:"<<fnvalue<<std::endl;
#else
    std::stringstream stst; std::string st="";
    for (int i=0;i<int(xv.size());i++) {  stst<<xv[i]; st+=stst.str()+" "; stst.str("");}
    st+=" f: "; stst<<fnvalue; st+=stst.str(); stst.str(""); Rprintf("%s\n",st.c_str());
#endif
    }
  } else resu=test->lambdaEst;
  return(wrap(resu));
}

//[[Rcpp::export]]
List Krig_coef_Wrapper(SEXP aA, SEXP lambdaP) {
  Rcpp::NumericVector xaA(aA);
  double lambda =as<double>(lambdaP);
  std::vector<covTypedef> CovTheta(0);
  for (Rcpp::NumericVector::iterator ii=xaA.begin();ii!=xaA.end()-1;ii++) CovTheta.push_back(*ii);
  //for (int ii=0;ii<xaA.size()-1;ii++) CovTheta.push_back(xaA(ii));
  test->Krig_engine_default<covTypedef,internalTypeEigen>(CovTheta, xaA(xaA.size()-1));
  test->Krig_coef<covTypedef,internalTypeEigen>(lambda);
  CKrigptrTable.push_back(test);
  // the position (C-style) of the CSmooth pointer in the table is passed back to R
  return(List::create(Named("u")=test->u,
                      Named("c")=test->coefs_random,Named("d")=test->coefs_fixed,
                      Named("D")=test->D_invEigVals,
                      Named("CKrigidx")=int(CKrigptrTable.size()-1)));
}

//[[Rcpp::export]]
int getFnEvalCount() { return(fnevalcounter); }

//[[Rcpp::export]]
SEXP CcovFocal( // ici la def
      SEXP focal,
      SEXP CKrigidxP) {
  int CKrigidx= as<int>(CKrigidxP);
  if (CKrigidx>=int(CKrigptrTable.size()) || CKrigidx<0) {
    {
      std::stringstream stst;stst<<"(!) Ccovfocal called with index out of allowed range, which is 0 -- "<<CKrigptrTable.size()-1<<std::endl;
#ifdef NO_R_CONSOLE
      std::cout<<stst.str();
#else
      REprintf(stst.str().c_str());
#endif
      throw Rcpp::exception("Ccovfocal called with index out of allowed range");
      //Rf_error("Ccovfocal called with index out of allowed range");
    }
  } else {
    NumericVector xfocal(focal);
    CSmooth* localCKrigptr=CKrigptrTable[CKrigidx];
    std::vector<ioType> focal_C(fittedparamnbr,0); // rmind that fittedparamnbr is a global variable that has been filled before;
    for (int ii=0;ii<fittedparamnbr;ii++) { focal_C[ii]=xfocal(ii); } // FR->FR avoid this copy if fillaxialFocal def is modified ?
    localCKrigptr->fillaxialFocal(focal_C);
    localCKrigptr->filleuclFocal(); //squared scaled distances
    localCKrigptr->fillcovFocal<covTypedef>();
    return(wrap(localCKrigptr->covFocal));
  }
}
