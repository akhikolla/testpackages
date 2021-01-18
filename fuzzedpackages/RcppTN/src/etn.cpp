/// Vectorized Truncated Normal RNG

///
/// This code wraps the etn1() C++ function to calculate the
/// expectation of the Truncated Normal distributions in a vectorized
/// fashion with Rcpp inputs and no output
///

# include <Rcpp.h>
# include "etn1.h"
# include "check1.h"

void etn(Rcpp::NumericVector &Mean, ///< vector of means
         Rcpp::NumericVector &Sd, ///< vector of standard deviations
         Rcpp::NumericVector &Low,	///< vector of lower bounds
         Rcpp::NumericVector &High,	///< vector of upper bounds
         Rcpp::NumericVector &Exps
         ) {

  // Namespace
  using namespace Rcpp ;
  //

  // Init Storage
  NumericVector::iterator itM = Mean.begin() ;
  NumericVector::iterator itS = Sd.begin() ;
  NumericVector::iterator itL = Low.begin() ;
  NumericVector::iterator itH = High.begin() ;
  NumericVector::iterator itE = Exps.begin() ;
  //

  // Expectations for TNs
  while (itE != Exps.end()) {
    if (check1(*itM, *itS, *itL, *itH)) {
      *itE = etn1(*itM, *itS, *itL, *itH);
    } else {
      *itE = NA_REAL ;
    }
      itE++ ;
      itM++ ;
      itS++ ;
      itL++ ;
      itH++ ;
  }
  //
}
