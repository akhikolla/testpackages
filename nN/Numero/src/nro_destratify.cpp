/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
RcppExport SEXP
nro_destratify(SEXP data_R, SEXP strata_R) {
  mdreal rlnan = medusa::rnan();
 
  /* Extract matrix data. */
  NumericMatrix mtx(data_R);
  mdsize nrows = mtx.nrow();
  mdsize ncols = mtx.ncol();
  
  /* Check inputs. */
  vector<mdsize> strata = nro::vector2sizes(strata_R);
  if(strata.size() != nrows)
    return CharacterVector("Incompatible inputs.");
  
  /* Process columns. */
  List output;
  for(mdsize j = 0; j < ncols; j++) {
    
    /* Check for non-finite values. */
    NumericVector values = mtx(_, j);
    LogicalVector flags = Rcpp::is_finite(values);
    
    /* Copy column vector. */
    vector<mdreal> array(nrows, rlnan);
    for(mdsize i = 0; i < nrows; i++) {
      if(!flags[i]) continue;
      array[i] = values[i];
    }
   
    /* Standardize distributions. */
    array = abacus::destratify(array, strata);
    
    /* Update results. */
    output.push_back(nro::reals2vector(array));
  }
  return output;
}
