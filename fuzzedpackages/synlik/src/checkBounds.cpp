#include "synlik.h"

/*
 *   NB cholFact_ is the _TRANSPOSE_ of the Cholesky factor
 */

SEXP checkBoundsCpp(SEXP theMean_, 
                    SEXP cholFact_, 
                    SEXP indexes_, 
                    SEXP upper_,
                    SEXP lower_,
                    SEXP output_)
{
  using namespace Rcpp;
  
  try{
    arma::colvec theMean = as<arma::colvec>(theMean_);
    arma::mat cholFact = as<arma::mat>(cholFact_);
    arma::mat output = as<arma::mat>(output_);
    
    NumericVector indexes = clone(as<NumericVector>(indexes_));
    NumericVector upper = as<NumericVector>(upper_);
    NumericVector lower = as<NumericVector>(lower_);
    
    int nsim = output.n_rows;
    int nPar = output.n_cols;
    int nChecks = indexes.size();
    indexes = indexes - 1.0;
    
    bool found = false;
    int jj;
    int iSimul = 0;
    
    // Check one row of output at the time
    while(iSimul < nsim)
    {
      found = false;
      jj = 0;
      
      // We move along each row as long as there are new elements to check (in indexes) or we find an element
      // that lays outside the boundaries.
      while( jj < nChecks && !found)
      {
        if( output(iSimul, indexes[jj]) > upper[jj] || output(iSimul, indexes[jj]) < lower[jj])
        {
          found = true;
        }
        jj++;
      }
      
      // If found == true the vector lays outside the boundaries and we simulate it again
      if(found) 
      {
        output.row(iSimul) = arma::trans( theMean + cholFact * as<arma::colvec>(rnorm(nPar)) );
      } else {
        iSimul++;
      }
      
    }
    return wrap(output);
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);
}

