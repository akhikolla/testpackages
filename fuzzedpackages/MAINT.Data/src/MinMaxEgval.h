#ifndef _Minmaxegval_h
#define _Minmaxegval_h

#define ARMA_DONT_PRINT_ERROR
#define ARMA_DONT_PRINT_WARNING
#include "RcppArmadillo.h"

bool MinMaxEgval(const arma::mat& Sigma, const int Cf, double& lnviol,
  double& minmppegv, double& maxmpegv, double& minlregv, double& maxlregv,  
  const double eps, const double minlndet, const double maxlnk2);

#endif

