/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 * Pearson's correlation coefficient.
 */
pair<mdreal, mdsize>
abacus::correlation(const vector<mdreal>& x, const vector<mdreal>& y) {
  mdsize i;
  mdsize nelem = x.size();
  mdreal rlnan = medusa::rnan();
  if(nelem != y.size())
    panic("Incompatible inputs.\n", __FILE__, __LINE__);

  /* Summary statistics. */
  mdsize nsum = 0;
  long double sx = 0.0;
  long double sy = 0.0;
  long double sxx = 0.0;
  long double syy = 0.0;
  long double sxy = 0.0;
  for(i = 0; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    if(y[i] == rlnan) continue;
    long double xval = x[i];
    long double yval = y[i];
    sx += xval;
    sy += yval;
    sxx += xval*xval;
    sxy += xval*yval;
    syy += yval*yval;
    nsum++;
  }

  /* Calculate coefficient. */
  pair<mdreal, mdsize> res(rlnan, nsum);
  long double vx = sqrt(nsum*sxx - sx*sx);
  long double vy = sqrt(nsum*syy - sy*sy);
  if(vx < 1e-20) return res;
  if(vy < 1e-20) return res;
  res.first = (nsum*sxy - sx*sy)/vx/vy;
  return res;
}
