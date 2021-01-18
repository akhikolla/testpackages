/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

#define NBURN 6
#define ALPHA 0.5

/*
 * Check the stability for a sequence of numbers.
 */
bool
abacus::convergence(const vector<mdreal>& history, const mdreal limit) {
  mdreal rlnan = medusa::rnan();

  /* Check the last three points. */
  mdsize nhist = history.size();
  if(nhist < 3) return false;
  mdreal last3 = history[nhist-3];
  mdreal last2 = history[nhist-2];
  mdreal last1 = history[nhist-1];
  if((last3 == last2) && (last3 == last1)) return true;
  if(nhist < NBURN) return false;
  
  /* Calculate differences. */
  mdsize nvalid = 0;
  mdreal prev = rlnan;
  long double asum = 0.0;
  long double dsum = 0.0;
  for(mdsize i = ALPHA*nhist; i < nhist; i++) {
    if(history[i] == rlnan) continue;
    mdreal x = prev;
    mdreal y = history[i];
    prev = history[i];
    if(x == rlnan) continue;
    if(y == rlnan) continue;
    asum += fabs(y - x);
    dsum += (y - x);
    nvalid++;
  }

  /* Convergence score. */
  if(nvalid < 1) return false;
  mdreal delta = fabs(dsum)/(asum + 1e-9);
  delta *= (nvalid + NBURN)/(nvalid + 1.0);
  return (delta < limit);
}
