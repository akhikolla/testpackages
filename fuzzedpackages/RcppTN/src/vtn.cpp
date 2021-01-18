/// Vectorized Truncated Normal RNG

///
/// This code wraps the vtn1() C++ function to calculate the
/// variance of the Truncated Normal distributions in a vectorized
/// fashion with Rcpp inputs and no output
///

# include <Rcpp.h>
# include "vtn1.h"
# include "check1.h"

void vtn(Rcpp::NumericVector &Mean, ///< vector of means
         Rcpp::NumericVector &Sd, ///< vector of standard deviations
         Rcpp::NumericVector &Low,	///< vector of lower bounds
         Rcpp::NumericVector &High,	///< vector of upper bounds
         Rcpp::NumericVector &Vars
         ) {

  // Namespace
  using namespace Rcpp ;
  //

  // Init Storage
  NumericVector::iterator itM = Mean.begin() ;
  NumericVector::iterator itS = Sd.begin() ;
  NumericVector::iterator itL = Low.begin() ;
  NumericVector::iterator itH = High.begin() ;
  NumericVector::iterator itV = Vars.begin() ;
  //

  // Expectations for TNs
  while (itV != Vars.end()) {
    if (check1(*itM, *itS, *itL, *itH)) {
      *itV = vtn1(*itM, *itS, *itL, *itH);
    } else {
      *itV = NA_REAL ;
    }
      itV++ ;
      itM++ ;
      itS++ ;
      itL++ ;
      itH++ ;
  }
  //
}
