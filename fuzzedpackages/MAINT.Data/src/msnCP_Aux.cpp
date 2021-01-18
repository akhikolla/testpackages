//#include "Rcpp.h"
#include "RcppArmadillo.h"
#include "msnCP_Aux.h"

using namespace Rcpp ;

const double LN2 = log(2.);
const double MINPHIARG = -50.;
// minimum value allowed for the argument of a cumulative distribution of a standardized random variable 

double zeta(const int k,const double x)  // k integer in (0,2)
{
  double x2,z1;
  static NumericVector x_ASNV(1),logphix_ASNV(1),logPhix_ASNV(1); 
  switch (k)  {
    case 0:
      x_ASNV(0) = x; 
      logPhix_ASNV = pnorm(x_ASNV,0.,1.,true,true); 
      return	logPhix_ASNV(0) + LN2;
    case 1:
      if (x>MINPHIARG) {
        x_ASNV(0) = x; 
        logphix_ASNV = dnorm(x_ASNV,0.,1.,true);
        logPhix_ASNV = pnorm(x_ASNV,0.,1.,true,true); 
        return exp(logphix_ASNV(0)-logPhix_ASNV(0));
      } else {
        x2 = x*x;
        return  -x / (	1 -1/(x2+2) +1/((x2+2)*(x2+4)) 
                          -5/((x2+2)*(x2+4)*(x2+6))
                          +9/((x2+2)*(x2+4)*(x2+6)*(x2+8)) 
                          -129/((x2+2)*(x2+4)*(x2+6)*(x2+8)*(x2+10))
                     );
      }
    case 2:
      z1 = zeta(1,x);
      return -z1*(x+z1);
      default: 
      return (double) NAN; 
  }
}


