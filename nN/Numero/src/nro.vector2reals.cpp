/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
vector<mdreal>
nro::vector2reals(const SEXP& data) {  
  mdreal rlnan = medusa::rnan();
  vector<mdreal> array;
  NumericVector values(data);
  LogicalVector flags = Rcpp::is_finite(values);
  mdsize nelem = values.size();
  for(mdsize i = 0; i < nelem; i++) {
    if(flags[i]) array.push_back(values[i]);
    else array.push_back(rlnan);
  }
  return array;
}
