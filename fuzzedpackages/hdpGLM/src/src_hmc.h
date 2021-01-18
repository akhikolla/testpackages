//=================================
// include guard
// -------------
#ifndef __HMC_H_INCLUDED__     // if .h hasn't been included yet...
#define __HMC_H_INCLUDED__     //   #define this so the compiler knows it has been included
//=================================



// forward declared dependencies
// ----------------------------
// class <class-to-include>

// Include dependencies
// --------------------
// #include <RcppArmadillo.h>

arma::colvec hmc_update(arma::colvec theta_t, double epsilon, int L, Rcpp::List fix);

//=================================
#endif
//=================================
