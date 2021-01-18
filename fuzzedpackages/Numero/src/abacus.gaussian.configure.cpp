/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
bool
Gaussian::configure(const vector<mdreal>& x, const vector<mdreal>& w) {
  mdreal rlnan = medusa::rnan();
  mdsize nelem = x.size();
  if(nelem < 2) return false;

  /* Sort data. */
  this->values = x;
  this->weights = w;
  vector<mdsize> sorted = sortreal(values, 1);
  for(mdsize k = 0; k < nelem; k++)
    weights[k] = w[sorted[k]];
  
  /* Determine location and scale. */
  this->offset = values[0];
  this->center = statistic(values, weights, "center");
  this->scale = (center - offset + 1e-16);

  /* Standardize data. */
  for(mdsize i = 0; i < nelem; i++)
    values[i] = (values[i] - offset)/scale;

  /* Calculate sum of weights. */
  long double wsum = 0.0;
  for(mdsize i = 0; i < nelem; i++)
    wsum += weights[i];

  /* Reference z-scores. */
  for(mdreal t = -10.0; t <= 10.0; t += 0.2)
    (this->zscores).push_back(t);

  /* Reference quantiles from normal distribution. */
  vector<mdreal> qref;
  mdreal qmin = weights[0]/wsum;
  for(mdsize i = 0; i < zscores.size(); i++) {
    mdreal z = zscores[i]/sqrt(2.0);
    mdreal q = 0.5*erfc(-z);
    if(q < qmin) q = rlnan;
    if(q >= 1.0) q = rlnan;
    qref.push_back(q);
  }

  /* Find evaluation elements. */
  long double wcount = 0.0;
  (this->qloci).resize(qref.size(), nelem);
  for(mdsize i = 0; i < (nelem - 1); i++) {
    mdreal a = wcount/wsum; wcount += weights[i];
    mdreal b = wcount/wsum;
    for(mdsize k = 0; k < qref.size(); k++) {
      mdreal q = qref[k];
      if((q == rlnan) || (a > q) || (b < q)) continue;
      qloci[k] = i;
    }
  }

  /* Exclude unusable elements. */
  mdsize nloci = 0;
  for(mdsize k = 0; k < zscores.size(); k++) {
    if(qloci[k] >= nelem) continue;
    zscores[nloci] = zscores[k];
    qloci[nloci] = qloci[k];
    nloci++;
  }

  /* Trim vectors. */
  (this->zscores).resize(nloci);
  (this->qloci).resize(nloci);
  
  /* Default approximation. */
  return (this->optimize("linear") >= 0.0);
}
