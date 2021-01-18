/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
NumericMatrix
nro::reals2matrix(const vector<vector<mdreal> >& vectors) {
  mdreal rlnan = medusa::rnan();
  mdsize nrows = vectors.size();
  if(nrows == 0) {
    NumericMatrix matrix(0, 0);
    return matrix;
  }
  NumericMatrix matrix(nrows, vectors[0].size());
  for(mdsize i = 0; i < vectors.size(); i++) {
    const vector<mdreal>& array = vectors[i];
    for(mdsize j = 0; j < array.size(); j++) {
      if(array[j] == rlnan)
        matrix(i, j) = NA_REAL;
      else
        matrix(i, j) = array[j];
    }
  }
  return matrix;
}
