// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

/// Vectorized Truncated Normal RNG

///
/// This code wraps the enttn1() C++ function to calculate the
/// entropy of the Truncated Normal distributions in a vectorized
/// fashion with Rcpp inputs and no output
///

# include <Rcpp.h>
# include "enttn1.h"

void enttn(Rcpp::NumericVector &Mean, ///< vector of means
           Rcpp::NumericVector &Sd, ///< vector of standard deviations
           Rcpp::NumericVector &Low,	///< vector of lower bounds
           Rcpp::NumericVector &High,	///< vector of upper bounds
           Rcpp::NumericVector &Ents
           ) {

    // Namespace
    using namespace Rcpp ;
    //

    // Init Storage
    NumericVector::iterator itM = Mean.begin() ;
    NumericVector::iterator itS = Sd.begin() ;
    NumericVector::iterator itL = Low.begin() ;
    NumericVector::iterator itH = High.begin() ;
    NumericVector::iterator itE = Ents.begin() ;
    //

    // Expectations for TNs
    while (itE != Ents.end()) {
        // if (check1(*itM, *itS, *itL, *itH)) {
        //   *itE = etn1(*itM, *itS, *itL, *itH);
        // } else {
        //   *itE = NA_REAL ;
        // }
        *itE = enttn1(*itM, *itS, *itL, *itH) ;
        itE++ ;
        itM++ ;
        itS++ ;
        itL++ ;
        itH++ ;
    }
    //
}
