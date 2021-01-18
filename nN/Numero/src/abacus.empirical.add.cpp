/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
bool
Empirical::add(const mdreal x, const mdreal w) {
  EmpiricalBuffer* p = (EmpiricalBuffer*)buffer;

  /* Check inputs. */
  mdreal rlnan = medusa::rnan();
  if(x == rlnan) return false;
  if(w == rlnan) return false;
  if(w <= 0.0) return false;

  /* Add value. */
  unordered_map<mdreal, mdreal>& data = p->data;
  if(data.count(x) < 1) data[x] = w;
  else data[x] += w;
  p->ndata += 1;

  /* Reset dependent data structures. */
  p->approx = Approximation();
  (p->valsorted).clear();
  (p->wsorted).clear();
  return true;
}
