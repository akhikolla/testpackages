/// R-accessible Vectorized Variance Calculation for Truncated Normals

///
/// This code is a wrapper for the vtn() C++ function to be called
/// from R.
///

# include <Rcpp.h>
# include "vtn.h"

RcppExport SEXP vtnRcpp(const SEXP mean, ///< vector of length K
					 ///containing the mean of the
					 ///distribution for each
					 ///calculation [real number]

			const SEXP sd, ///< vector of length K
				       ///containing the standard
				       ///deviation for each
				       ///calculation [strictly
				       ///positive real number]
			const SEXP low, ///< vector of length K
					///containing the lower bound
					///for each calculation [real
					///number or '-Inf']
			const SEXP high ///< vector of length K
					///containing the upper bound
					///for each calculation [real
					///number or 'Inf']
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
  NumericVector Vars(Mean.size(), 0.0) ;
  vtn(Mean, Sd, Low, High, Vars) ;
  return(Vars);
  //
}
