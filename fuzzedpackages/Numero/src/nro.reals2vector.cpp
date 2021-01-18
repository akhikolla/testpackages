/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
NumericVector
nro::reals2vector(const vector<mdreal>& values) {
  mdreal rlnan = medusa::rnan();
  mdsize nelem = values.size();
  NumericVector array(nelem);
  for(mdsize i = 0; i < nelem; i++) {
    if(values[i] == rlnan) array[i] = NA_REAL;
    else array[i] = values[i];
  }
  return array;
}
