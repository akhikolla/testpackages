/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "nro.h"

/*
 *
 */
vector<mdsize>
nro::vector2sizes(const SEXP& data) {  
  mdsize sznan = medusa::snan();
  vector<mdsize> array;
  NumericVector values(data);
  LogicalVector flags = Rcpp::is_finite(values);
  mdsize nelem = values.size();
  for(mdsize i = 0; i < nelem; i++) {
    if(flags[i]) array.push_back((mdsize)(values[ i ] + 0.5));
    else array.push_back(sznan);
  }
  return array;
}
