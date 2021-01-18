#ifndef H_KRIGMAIN
#define H_KRIGMAIN

#include <cstdlib>
#include <string>
#include <sstream>
#include <limits>
#include <cmath> //for std::abs
#include "Rmath.h" // bessel_k (and apparently required before loading Rcpp.h...)
#include "Bessel_nr.h" // bessk
#include <Rcpp.h> // a charger apres les headers contenant des templated functions ???
//                  (sinon par exemple il ne comprend plus error())

namespace NS_GG {
 extern int a;
 extern covTypedef b;
};


bool newCSmooth( //header declaration
    Rcpp::NumericMatrix xy,
    int nrowxy, //with replicates
    int ncolxy,
    int nuniquerows, // required to tell the allocated size of c, d, D, u in *R*
    //    double *maxSmoothness,
    //    double *c,
    //    double *d,
    //    double *D,
    //    double *u,
    //    double *lambda,
    double GCV,
    //    double *covfnparam, //p+1 form
    //int fneval,
    int optimiseBool,
    int verbosity
    //    double *initcovfnparam //p+1 form  // new arg 2015/10/24
);



#endif

