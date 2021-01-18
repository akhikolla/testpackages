/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 * Interpolate Y at T. Any T values outside X range are not estimated
 * (not-a-number). If X contains duplicate values, the choice of the
 * corresponding Y is arbitrary.
 */
vector<mdreal>
abacus::interpolate(const vector<mdreal>& x0, const vector<mdreal>& y0,
		    const vector<mdreal>& t) {
  mdsize sznan = medusa::snan();
  mdreal rlnan = medusa::rnan();
  mdsize nelem = x0.size();

  /* Check inputs. */
  if(nelem < 1) panic("No data.", __FILE__, __LINE__);
  if(y0.size() != nelem)
    panic("Incompatible inputs.", __FILE__, __LINE__);

  /* Check if already sorted. */
  bool sortflag = false;
  vector<mdreal> x = x0;
  vector<mdreal> y = y0;
  for(mdsize i = 1; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    if(y[i] == rlnan) continue;
    if(x[i] >= x[i-1]) continue;
    sortflag = true;
    break; 
  }

  /* Sort coordinates. */
  if(sortflag) {
    vector<mdsize> sorted = sortreal(x, 1);
    for(mdsize k = 0; k < nelem; k++)
      y[k] = y0[sorted[k]];
  }
  
  /* Remove duplicate points. */
  mdsize nvalid = 0;
  for(mdsize i = 0; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    if(y[i] == rlnan) continue;
    if(i > 0) if(x[i] == x[i-1]) continue;
    x[nvalid] = x[i];
    y[nvalid] = y[i];
    nvalid++;
  }

  /* Trim arrays. */
  x.resize(nvalid);
  y.resize(nvalid);
  if(nvalid < 2) panic("Not enough data.", __FILE__, __LINE__);

  /* Prepare output array. */
  mdsize nout = t.size();
  if(nout < 1) return t;
  vector<mdreal> z(nout, rlnan);

  /* Interpolate at new sampling points. */
  mdreal xmin = x[0];
  mdreal xmax = x[nvalid-1];
  for(mdsize i = 0; i < t.size(); i++) {
    mdreal key = t[i];
    if(key == rlnan) continue;
    if(key < xmin) continue;
    if(key > xmax) continue;

    /* Interpolate between nearest points. */
    Site slot = binsearch(x, key);
    if(slot.bounds.first == sznan) continue;
    if(slot.bounds.second == sznan) continue;
    z[i] = (slot.weights.first)*(y[slot.bounds.first])
      + (slot.weights.second)*(y[slot.bounds.second]);
  }
  return z;
}
