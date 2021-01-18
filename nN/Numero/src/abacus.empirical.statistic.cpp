/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Empirical::statistic(const std::string& name) const {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;
  vector<mdreal> values;
  vector<mdreal> weights;
  p->contents(values, weights);
  return abacus::statistic(values, weights, name);
}
