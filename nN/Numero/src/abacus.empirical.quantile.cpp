/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Empirical::quantile(const mdreal q) const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  vector<mdreal> values;
  vector<mdreal> weights;
  p->contents(values, weights);
  return abacus::quantile(values, weights, q);
}
