/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdsize
BaseGaussian::transform(vector<mdreal>& x) const {
  mdreal rlnan = medusa::rnan();

  /* No model. */
  if(sigma == rlnan) {
    for(mdsize i = 0; i < x.size(); i++)
      x[i] = rlnan;
    return 0;
  }

  /* Normalize data. */
  for(mdsize i = 0; i < x.size(); i++) {
    if(x[i] == rlnan) panic("Unusable value.", __FILE__, __LINE__);
    x[i] = (x[i] - offset)/scale;
  }

  /* Apply transformation. */
  this->apply(x, factor);

  /* Estimate Z-scores. */
  for(mdsize i = 0; i < x.size(); i++)
    x[i] = (x[i] - mu)/sigma;
  return x.size();
}
