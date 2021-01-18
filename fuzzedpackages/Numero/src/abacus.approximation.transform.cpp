/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Approximation::transform(const mdreal x) const {
  mdreal rlnan = medusa::rnan();
  if(x == rlnan) return rlnan;
  if(mode == rlnan) return rlnan;

  /* Two-sided transform. */
  vector<mdreal> vals(1);
  mdreal zP; vals[0] = x; positive.transform(vals); zP = vals[0];
  mdreal zN; vals[0] = -x; negative.transform(vals); zN = -vals[0];
  if(zP == rlnan) return rlnan;
  if(zN == rlnan) return rlnan;

  /* Combination weight coefficients. */
  mdreal rP = 0.5;
  mdreal rN = 0.5;
  if(x > mode) {
    rN = 0.5/(1.0 + fabs(zP));
    rP = (1.0 - rN);
  }
  if(x < mode) {
    rP = 0.5/(1.0 + fabs(zN));
    rN = (1.0 - rP);
  }

  /* Interpolate final score. This fixes gaps around the mode. */
  return (rN*zN + rP*zP);
}
