/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 * The input 'f' should be within ]0, 1[.
 */
void
BaseGaussian::apply(vector<mdreal>& x, const mdreal f) const {
  if(f < 0.0) panic("Unusable input.", __FILE__, __LINE__);
  if(f > 1.0) panic("Unusable input.", __FILE__, __LINE__);

  /* Non-linear transforms. */
  if(method == "exp") {
    mdreal g = 7*f*f;
    for(vector<mdreal>::iterator it = x.begin(); it != x.end(); it++)
      *it = exp(g*(*it));
    return;
  }
  if(method == "log") {
    mdreal g = 8*(f - 0.5); g = exp(g*fabs(g));
    for(vector<mdreal>::iterator it = x.begin(); it != x.end(); it++) {
      mdreal value = (*it/g + 1.0);
      if(value <= 1e-20) value = 1e-20;
      *it = log(value);
    }
    return;
  }

  /* No transform. */
  if(method == "linear") return;
  panic("Unknown method.", __FILE__, __LINE__);
}
