/// R-accessible Vectorized Truncated Normal RNG

///
/// This code is a wrapper for the rtn() C++ function to be called
/// from R.
///

# include <Rcpp.h>
# include "rtn.h"

RcppExport SEXP rtnRcpp(const SEXP mean, ///< vector of length K
					 ///containing the mean of the
					 ///distribution for each draw
					 ///[real number]

			const SEXP sd, ///< vector of length K
				       ///containing the standard
				       ///deviation for each draw
				       ///[strictly positive real
				       ///number]
			const SEXP low, ///< vector of length K
					///containing the lower bound
					///for each draw [real number
					///or '-Inf']
			const SEXP high ///< vector of length K
					///containing the upper bound
					///for each draw [real number
					///or 'Inf']
			) {

  // Namespace
  using namespace Rcpp ;
  //

  // Conversion of Inputs
  NumericVector Mean(mean) ;
  NumericVector Sd(sd) ;
  NumericVector Low(low) ;
  NumericVector High(high) ;
  //

  // Init, Populate, and Return
  NumericVector Draws(Mean.size(), 0.0) ;
  {
    RNGScope scope ;
    rtn(Mean, Sd, Low, High, Draws) ;
  }
  return Draws ;
  //
}
