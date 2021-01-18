/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
Approximation::fit(const vector<mdreal>& x, const vector<mdreal>& w) {
  mdreal rlnan = medusa::rnan();

  /* Check if anything to do. */
  if(mode != rlnan) return;
  
  /* Find the center of the distribution. */
  this->mode = statistic(x, w, "mode");

  /* Create a mirror version for decreasing tail. */
  vector<mdreal> xN = x;
  for(mdsize i = 0; i < xN.size(); i++)
    if(xN[i] != rlnan) xN[i] *= -1;

  /* Simple approximation. */
  (this->positive).configure(x, w);
  (this->negative).configure(xN, w);

  /* Determine the best model for increasing tail. */
  Gaussian gauss = this->positive;
  if(gauss.optimize("exp") > positive.quality()) this->positive = gauss;
  if(gauss.optimize("log") > positive.quality()) this->positive = gauss;

  /* Determine the best model for decreasing tail. */
  gauss = this->negative;
  if(gauss.optimize("exp") > negative.quality()) this->negative = gauss;
  if(gauss.optimize("log") > negative.quality()) this->negative = gauss;
}
