
#include "synlik.h"

/* Takes an input matrix X and checks if there are Nas or NaN row by row
* if there some rows that contain at least 1 NA or NaN
* it returns a list of 
* 1) Total number of contaminated rows
* 2) Indexes of the contaminated rows
* 3) a new "clean" matrix where the bad rows have been removed
*/

SEXP cleanStats(SEXP inMat)
{
    using namespace Rcpp;
    
    try{
    NumericMatrix X = as<NumericMatrix>(inMat);
    
    int nRows = X.nrow();
    int nCols = X.ncol();
    
    //if(nCols > nRows) Rcout << "cleanStats:  nCols > nRows mean more statistics then simulations?!" << std::endl;
    
    IntegerVector banned;
    
    // Looking row by row...      
    for(int iRow = 0; iRow < nRows; iRow++)
      for(int iCol = 0; iCol < nCols; iCol++)
        if(R_IsNA(X(iRow, iCol)) || R_IsNaN(X(iRow, iCol)))
        {
          banned.push_back(iRow); 
          iCol = nCols;            // We stop checking this row is we find an NA or NaN
        }
    
    int nBanned = banned.size();
    if(nBanned == nRows) stop("All the vectors of summaries statistics contain NAs or NaNs");
    
    /* IF( there are some NaN or NA ) we fill up a new "clean" matrix row by row.
    * ELSE we return a list zeros
    * The following code is quite messy! 
      */
      if(nBanned > 0)
      {
        int remRows = nRows - nBanned;
        NumericMatrix cleanX(remRows, nCols);
        
        int banIndex = 0;
        int outRow = 0;
        for(int inRow = 0; inRow < remRows; outRow++)
        {
          if(banIndex >= nBanned || outRow != banned[banIndex])
          {
            cleanX(inRow, _) = X(outRow, _);
            inRow++;
          } else {
            banIndex++; 
          } 
        }
        
        return List::create(_["nBanned"] = nBanned,
                            _["banned"] = (banned + 1), // +1 because we return R indexes
                            _["cleanX"] = cleanX  );
      } else {
        return List::create(_["nBanned"] = 0,
                            _["banned"] = 0, 
                            _["cleanX"] = 0);
      }
    } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);  
}

  
