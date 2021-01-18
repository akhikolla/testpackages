/*====================================================================================
 Pijavski one-dimentional global minimization algorithm for non-smooth multiextremal functions

Parameters:
 SEXP Pijavski(SEXP o1, SEXP o2, SEXP fn, SEXP Lips, SEXP x1, SEXP x2, SEXP iter, SEXP prec, SEXP rho)

 inputs: fn - pointer to the objectibe function
 Lips - Lipschitz constant
 a,b  - the interval on which to optimize
 iter - limit on the number of iterations (on return, number of iterations made)
 prec - desired precision (in terms of the value of the objective the difference between obj and global lower bound on fn)
   if negative on return, then the Lipschitz constant is too small
  
 env - environment variable passed from R, it should be defined as new.env(list(fn = myfunction))
       where myfunction is the objective function


returns the list
  (x,val,precision, maxiter)
  x - the global minimum
  val - value at the minimum
  precision is precision reached (if negative, then Lipschitz constant is too small, repeat with a larger Lips)
  maxiter is the number of iterations made
 
 Freeware. Written by Gleb Beliakov, 2004
 gleb@dekain.edu.au 
 Packaged for R by Gita Das, 2016 

=======================================================================================*/
#include "heap.h"

#include <cstdlib>
#include <vector>

#include <Rcpp.h>
#include <Rdefines.h>
#include <stdint.h>

using namespace Rcpp;


//typedef void ( *USER_FUNCTION_P)(double *, double *);

//The folowing lines are used to disable optimizations since GCC 4.4
// #pragma GCC push_options
// #pragma GCC optimize("O0")
// this generates warnings for GCC ver 4.6.3




#ifdef _MSC_VER
// if the compiler does not recognise this type, change it to another int type 8 bytes long
// like long long int
typedef  __int64 ULINT; //this type myst be 8 bytes long
#else
// #pragma GCC diagnostic ignored "-Wno-long-long"
// typedef uint_least64_t ULINTH;
typedef unsigned long long int ULINT; //this type myst be 8 bytes long
#endif


KEY_TYPE ftP; // global

#define PACKFLOAT(f) ( (f  & 0xFFFFFFFF))
//the following pragma is used to avoid compiler warnings
// "dereferencing type-punned pointer will break strict-aliasing rules" as broken in the two macros
// f2ulint and ulint2f (both pointer int and float refering to same memory locations )
//#pragma GCC diagnostic ignored "-Wstrict-aliasing"

/*  GB sacrifice space over security/compatibility
#define f2ulint(a)  (*(unsigned int*) &(a))
#define ulint2f(a) (*(float*) &(a))
unsigned int UL1P;

KEY_TYPE _getKey(node_t theNode)
{
	UL1P=theNode.data>>32;
	ftP = ulint2f(UL1P);
	// use the first 4 bytes
	return -ftP;
}
*/

INDEX_TYPE _getIndex(node_t theNode) {	return 0; }
void _setIndex(node_t theNode, INDEX_TYPE I) {}

KEY_TYPE _getKey(node_t theNode)
{
	return -(theNode.data.val);
}


ULINT ULP;
unsigned int UI;
DATA_TYPE DAT;
/*
INDEX_TYPE Merge1 (bheap_t *h, float f, unsigned short i, unsigned short j)
{
	f=-f;
	ULP=f2ulint(f);//<<32;

#ifdef _MSC_VER
       ULP = (ULP << 32) & 0xFFFFFFFF00000000UL ;
#else
       ULP = (ULP << 32) & 0xFFFFFFFF00000000ULL ;
#endif
    UI = i;
	UI = ((UI  << 16) & 0xFFFF0000) | j;
	ULP = ULP + UI;
	return bh_insert(h, ULP);
}*/
INDEX_TYPE Merge1(bheap_t *h, float f, unsigned int  i, unsigned int j)
{
	f = -f;
	DAT.val = f;
	DAT.Idx1 = i;
	DAT.Idx2 = j;
	return bh_insert(h, DAT);
}

/*
void GetIndices(DATA_TYPE Node, unsigned int *i, unsigned int *j)
{
#ifdef _MSC_VER
       UI = Node & 0x00000000FFFFFFFFUL ;
#else
       UI = Node & 0x00000000FFFFFFFFULL ;
#endif
	   *j = UI & 0x0000FFFF;
	   *i = (UI >> 16) & 0x0000FFFF;
}
*/
void GetIndices(DATA_TYPE Node, unsigned int *i, unsigned int *j)
{
	*i = Node.Idx1;
	*j = Node.Idx2;
}

//#pragma GCC optimize( "", on )
//#pragma GCC pop_options
//#pragma GCC diagnostic warning "-Wstrict-aliasing"

void ComputeMin(double x1, double x2, double f1, double f2, double M, double* t, float* f)
{
	*t= 0.5*(x1+x2) + 0.5/M*(f1-f2);
	*f= 0.5*(f1+f2) + 0.5*M*(x1-x2);
}

//' Pijavski code in C++ is being called from R environment with multiple arguments
//' It then periodically calls a function in R
//' "//[[Rcpp::export]]" is needed to export the function to R environment
//' @param str input character vector
//' @exportPattern("^[[:alpha:]]+")
//' @useDynLib Pijavski
//' @importFrom Rcpp evalCpp
//[[Rcpp::export]]
SEXP Pijavski( SEXP fn, SEXP Lips, SEXP x1, SEXP x2, SEXP iter, SEXP prec, SEXP rho) {
  SEXP fval, env = rho;
  double *x0, *val, *Lip, *Xl, *Xu, *precision;
  double *maxiter;

  //typedef void ( *USER_FUNCTION_P)(double *, double *);
	double o1=0,o2=0;
 // x0 = NUMERIC_POINTER(o1);
 // val = NUMERIC_POINTER(o2);
  
  x0=&o1;
  val=&o2;
  Lip = NUMERIC_POINTER(Lips);
  Xl = NUMERIC_POINTER(x1);
  Xu = NUMERIC_POINTER(x2);
  maxiter = NUMERIC_POINTER(iter);
  precision = NUMERIC_POINTER(prec);

  //Rprintf("%s", "Within Pijavski-C\n");
  //Rprintf("*Lip: %f,*x0: %f, *val: %f,  *Xl: %f, *Xu: %f, *maxiter: %f, *precision: %f", *Lip,*x0, *val,  *Xl, *Xu,*maxiter, *precision);

  bheap_t* HeapP = bh_alloc();
  DATA_TYPE Node;
  std::vector<double> x, f;

 // if (*maxiter>=0xFFFD) *maxiter=0xFFFD; // cannot use short indices!

  int Iter=0;
  double Best=10e10;
  double CurrPrecision=10e10;

  double t,t1=0, t2=0,t3=0,f1=0,f2=0, f3=0, M=*Lip;
  float ff;

  unsigned int i,j,k;

  t1=*Xl;
  //F(&t1,&f1);
  Rcpp::Function func(fn);
  PROTECT(fval = Rf_eval(func(t1, f1), env)); //evaluates the function in R, fval is the output
  //Rprintf("\n%s", "Back to Pijavski-C: ");
  f1 = Rf_asReal(fval);
  if(Best>f1) { Best=f1; *x0=t1;}
  //Rprintf("\nBest:%f, *x0:%f\n", Best,*x0);
  x.push_back(t1);
  f.push_back(f1);
  UNPROTECT(1);
  
  t2=*Xu;
  //F(&t2,&f2);
  PROTECT(fval = Rf_eval(func(t2,f2), env));
  f2 = Rf_asReal(fval);
  UNPROTECT(1);
  if(Best>f2) { Best=f2; *x0=t2;}
  x.push_back(t2);
  f.push_back(f2);

  i=0; j=1;
  ComputeMin(t1,t2,f1,f2,M,&t,&ff);
  CurrPrecision = Best - ff;

  Merge1(HeapP, -ff, i,j);

  Iter=1;
  while(Iter < *maxiter && CurrPrecision > *precision ) {
    Iter++;
    Node=bh_delete_min(HeapP);
    GetIndices(Node, &i, &j);

    t1=x[i]; t2=x[j]; f1=f[i]; f2=f[j];

    ComputeMin( t1,t2,f1,f2,  M, &t3, &ff);

    if(ff>= *val) // cannot achieve required minimum
    {
      goto exitL;
    }

    //F(&t3,&f3);
    PROTECT(fval = Rf_eval(func(t3,f3), env));
    f3 = Rf_asReal(fval);

    x.push_back(t3);
    f.push_back(f3);
    k=Iter;

    if(Best>f3) { Best=f3; *x0=t3;}
    //Rprintf("\nIter:%d, Best:%f, f3:%f, t3:%f, *x0:%f\n", Iter,Best, f3, t3,*x0);
    CurrPrecision = Best - ff;

    // two new minima
    ComputeMin(t1,t3,f1,f3,M,&t,&ff);
    Merge1(HeapP, -ff, i,k);

//    ComputeMin(t3,t2,f3,f2,M,&t,&ff);  // always the same value ff
    Merge1(HeapP, -ff, k,j);
    UNPROTECT(1);
  }

exitL:
    *val=Best;

    bh_free(HeapP);
    x.clear(); f.clear();

    *precision=CurrPrecision;
    *maxiter=Iter;
    if(CurrPrecision<0) return (Rf_ScalarInteger(-1));  // too small Lip const

    return (Rcpp::List::create(*x0,*val,*precision, *maxiter) );

  //all C functions have SEXP inputs and returns a SEXP
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation. R code should start with "/***R" and end with "*/"
//

/*** R
optimize_funcR <- function(x,y){
  y <- x * x
  return(y)
}

output<-Pijavski(optimize_funcR, 5, -2.0, 1.0, 1000, 10^-3, new.env(list(fn = optimize_funcR)))
output

*/
