#ifndef __dePril_h__
#define __dePril_h__

#include <RcppArmadillo.h>

arma::vec getProbs_dePril(unsigned xnum, const Rcpp::List distPars, 
			  arma::vec extrapolPars, const std::string dist,
			  const unsigned& nsteps, double time,
			  bool extrap);

arma::vec getProbs_dePril(unsigned xnum, const Rcpp::List distPars, 
			  arma::vec extrapolPars, Rcpp::Function survR,
			  const unsigned& nsteps, double time, bool extrap);

#endif
