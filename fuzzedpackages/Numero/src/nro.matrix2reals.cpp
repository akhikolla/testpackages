/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 * Convert column-major R-matrix into C++ row vectors.
 */
vector<vector<mdreal> >
nro::matrix2reals(const SEXP& data, const mdreal jitter) {  
  mdreal rlnan = medusa::rnan();

  /* Extract matrix data. */
  NumericMatrix mtx(data);
  mdsize nrows = mtx.nrow();
  mdsize ncols = mtx.ncol();
  
  /* Convert to C++ row matrix. */
  vector<vector<mdreal> > rows(nrows);
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

    /* Add pseudo-random jitter. */
    if(jitter > 0.0) {
      mdreal sigma = abacus::statistic(array, "sd");
      if(sigma == rlnan) sigma = 0.0;
      for(mdsize k = 0; k < array.size(); k++) {
	if(array[k] == rlnan) continue;
	mdreal amp = 0.05*((k*77 + j*19)%39 - 19.0);
	array[k] += jitter*sigma*amp;
      }
    }

    /* Update rows. */
    for(mdsize i = 0; i < nrows; i++)
      rows[i].push_back(array[i]);
  }
  return rows;
}
