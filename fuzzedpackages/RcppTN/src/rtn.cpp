/// Vectorized Truncated Normal RNG

///
/// This code wraps the rtn1() C++ function to call the Truncated
/// Normal RNG multiple times in a vectorized fashion with Rcpp inputs
/// and no output
///

# include <Rcpp.h>
# include "rtn1.h"
# include "check1.h"

void rtn(Rcpp::NumericVector &Mean, ///< vector of means
         Rcpp::NumericVector &Sd, ///< vector of standard deviations
         Rcpp::NumericVector &Low,	///< vector of lower bounds
         Rcpp::NumericVector &High,	///< vector of upper bounds
         Rcpp::NumericVector &Draws
         ) {

  // Namespace
  using namespace Rcpp ;
  //

  // Init Storage
  NumericVector::iterator itM = Mean.begin() ;
  NumericVector::iterator itS = Sd.begin() ;
  NumericVector::iterator itL = Low.begin() ;
  NumericVector::iterator itH = High.begin() ;
  NumericVector::iterator itD = Draws.begin() ;
  //

  // Draw from TN
  while (itD != Draws.end()) {
    if (check1(*itM, *itS, *itL, *itH)) {
      *itD = rtn1(*itM, *itS, *itL, *itH);
    } else {
      *itD = NA_REAL ;
    }
      itD++ ;
      itM++ ;
      itS++ ;
      itL++ ;
      itH++ ;
  }
  //
}
