//
//  'includes' and global variable definitions for class Clmbr


#if !defined  CLMBR_G_H__		// prevents compiler from repeating this code in other files
#define  CLMBR_G_H__


#include <iomanip>				// for 'setw'
#include <algorithm>			// for 'min', 'max'
#include <cmath>				// for 'sqrt', 'log', 'cos' 
#include <ctime>


#include <R.h>
#include <R_ext/Applic.h>		// for 'Rdqags' and 'Rdqagi'
#include <R_ext/Lapack.h>


#include <Rcpp.h>

#include "tnt_vector.h"


#ifdef ENABLE_NLS
   #include <libintl.h>
   #define _(String) dgettext ("lm.br", String)
#else
   #define _(String) (String)
#endif


using TNT::Vector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::stop;
using Rcpp::Rcout;
using std::endl;


enum  MODEL { M1, M2, M3 };
enum  METHOD { GEO, GEO2, AF, AF2, MC, INIT };
const double zero_eq = ldexp( 1., -40 );
const double Inf = R_PosInf;
const double NaN = R_NaN;
const double NA = NA_REAL;
const double pi = M_PI;


#endif


