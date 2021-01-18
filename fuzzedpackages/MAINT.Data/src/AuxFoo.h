#ifndef _AuxFoo_h
#define _AuxFoo_h

#include "RcppArmadillo.h"

using namespace Rcpp ;
using namespace arma ;

void outerprod(const int p,const vec& v1,const vec& v2,mat& res);
void outerprod(const int p,const vec& v,mat& res);
// mat RestCov(const int q,NumericVector::iterator xpos,const int Config,const bool FixedArrays);
mat RestCov(const int q,NumericVector::iterator xpos,const int Config,const bool FixedArrays, const bool Srpar);

#endif
