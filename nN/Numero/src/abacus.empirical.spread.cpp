/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdsize
Empirical::spread() const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  return (p->data).size();
}

/*
 *
 */
mdsize
Empirical::spread(vector<mdreal>& values, vector<mdreal>& weights) const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  p->contents(values, weights);
  return values.size();
}
