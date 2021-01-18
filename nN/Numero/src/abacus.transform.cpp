/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

static vector<mdreal> rank_tf(const vector<mdreal>&, const bool);

/*
 *
 */
vector<mdreal>
abacus::transform(const vector<mdreal>& x, const string& name) {
  mdsize nelem = x.size();
  mdreal rlnan = medusa::rnan();
  vector<mdreal> y;

  /* Calculate t-statistics (sample z-scores). */
  if(name == "z") {
    mdreal mu = statistic(x, "mean");
    mdreal sigma = statistic(x, "sd");
    if(sigma == rlnan) return x;
    if(sigma < 1e-20) return vector<mdreal>(nelem, 0.0);
    y.resize(nelem, rlnan);
    for(mdsize i = 0; i < nelem; i++) {
      if(x[i] == rlnan) continue;
      y[i] = (x[i] - mu)/sigma;
    }
    return y;
  }

  /* Rank-based transformations. */
  if(name == "rank") return rank_tf(x, false);
  if(name == "uniform") return rank_tf(x, true);
  if(name == "balanced") {
    y = rank_tf(x, true);
    for(mdsize i = 0; i < y.size(); i++) {
      if(y[i] == rlnan) continue;
      y[i] = (2.0*y[i] - 1.0);
    }
    return y;
  }
  if(name == "tapered") {
    y = rank_tf(x, true);
    for(mdsize i = 0; i < y.size(); i++) {
      if(y[i] == rlnan) continue;
      mdreal t = (2.0*y[i] - 1.0);
      y[i] = (t + pow(t, 3.0) + pow(t, 5.0))/3;
    }
    return y;
  }

  /* Invalid argument. */
  panic("Unknown transformation.\n", __FILE__, __LINE__);
  return y;
}

/*
 *
 */
vector<mdreal>
rank_tf(const vector<mdreal>& x, const bool flag) {
  mdsize nelem = x.size();
  mdreal rlnan = medusa::rnan();

  /* Sort data vector. */
  vector<mdreal> y = x;
  vector<mdsize> mask = medusa::sortreal(y, 1);

  /* Convert to ranks. */
  mdsize nvalid = 0;
  vector<mdsize> r(nelem, nelem);
  for(mdsize i = 0; i < nelem; i++) {
    if(y[i] == rlnan) continue;
    r[i] = nvalid;
    nvalid++;
  }

  /* Check if enough data. */
  if(nvalid < 1) return x;
  if(nvalid < 2) {
    for(mdsize i = 0; i < nelem; i++) {
      y[i] = rlnan;
      if(x[i] != rlnan) y[i] = 0.5;
    }
    return y;
  }

  /* Add a sentinel element. */
  y.push_back(rlnan);

  /* Adjust for duplicates. */
  for(mdsize i = 0; i < nelem; i++) {
    if(y[i] == rlnan) continue;

    /* Count duplicates. */
    mdsize n;
    for(n = 0; (i + n) <= nelem; n++) 
      if(y[i+n] != y[i]) break;

    /* Singular value. */
    if(n < 2) {
      y[i] = r[i];
      continue;
    }

    /* Duplicates get the same rank. */
    mdreal mu = (0.5*r[i] + 0.5*r[i+n-1]);
    for(mdsize k = i; k < (i + n); k++)
      y[k] = mu;
    i += (n - 1);
  }
 
  /* Normalize ranks. */
  if(flag) {
    for(mdsize i = 0; i < nelem; i++) {
      if(y[i] == rlnan) continue;
      y[i] /= (nvalid - 1.0);
      if(y[i] < 0.0) y[i] = 0.0;
      if(y[i] > 1.0) y[i] = 1.0;
    }
  }

  /* Restore original order. */
  vector<mdreal> z(nelem);
  for(mdsize i = 0; i < nelem; i++)
    z[mask[i]] = y[i];
  return z;
}
