/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
abacus::quantile(const vector<mdreal>& data, const mdreal q) {
  mdsize nelem = data.size();
  mdreal rlnan = medusa::rnan();
  if(nelem < 1) return rlnan;  
  if(q < 0.0) return rlnan; 
  if(q > 1.0) return rlnan; 

  /* Count the number of available data. */
  mdsize n = 0;
  vector<mdreal> x = data;
  for(mdsize i = 0; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    x[n] = x[i];
    n++; 
  }

  /* Trim array. */
  x.resize(n);
  if(n < 1) return rlnan;  

  /* Sort array. */
  sort(x.begin(), x.end());

  /* Determine quantile indices. */
  mdsize a = (mdsize)(q*(n - 1));
  mdsize b = (a + 1);
  if(b >= n) b = (n - 1);
  if(x[a] == x[b]) return x[a];

  /* Interpolate between indices. */
  double d = (q*(n - 1) - a);
  return ((1.0 - d)*x[a] + d*x[b]);
}

/*
 *
 */
mdreal
abacus::quantile(const vector<mdreal>& data,
		 const vector<mdreal>& weights, 
		 const mdreal q) {
  mdsize nelem = data.size();
  mdreal rlnan = medusa::rnan();
  mdsize sznan = medusa::snan();
  if(nelem < 1) return rlnan;  
  if(q < 0.0) return rlnan; 
  if(q > 1.0) return rlnan; 

  /* Extreme values. */
  if(q == 0.0) {
    mdsize ind = extrema(data).first;
    if(ind == sznan) return data[ind];
    return rlnan;
  }
  if(q == 1.0) {
    mdsize ind = extrema(data).second;
    if(ind == sznan) return data[ind];
    return rlnan;
  }

  /* Check weights. */
  vector<mdreal> w = weights;
  if(w.size() != nelem)
    panic("Incompatible inputs.", __FILE__, __LINE__);
  if(statistic(w, "range") == 0.0)
    return abacus::quantile(data, q);

  /* Count the number of available data. */
  mdsize nvalid = 0;
  vector<mdreal> x = data;
  for(mdsize i = 0; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    if(w[i] == rlnan) continue;
    if(w[i] <= 0.0) continue;
    x[nvalid] = x[i];
    w[nvalid] = w[i];
    nvalid++;
  }

  /* Trim arrays. */
  x.resize(nvalid);
  w.resize(nvalid);
  if(nvalid < 1) return rlnan;  
  if(nvalid == 1) return x[0];

  /* Sort array and estimate cumulative density. */
  vector<mdreal> t((nvalid + 1), 0.0);
  vector<mdreal> f((nvalid + 1), 0.0);
  vector<mdsize> sorted = sortreal(x, 1);
  for(mdsize k = 0; k < nvalid; k++) {
    mdsize ind = sorted[k];
    f[k+1] = (f[k] + w[ind]);
    t[k+1] = (k + 1.0)/nvalid;
  }

  /* Normalize cumulative density. */
  for(mdsize i = 0; i <= nvalid; i++)
    f[i] /= f[nvalid];

  /* Interpolate cumulative density back to original size. */
  vector<mdreal> u(nvalid);
  for(mdsize i = 0; i < nvalid; i++)
    u[i] = i/(nvalid - 1.0);
  f = interpolate(t, f, u);

  /* Determine quantile indices. */
  Site slot = binsearch(f, q);
  mdsize a = slot.bounds.first;
  mdsize b = slot.bounds.second;
  if(a == sznan) a = b;
  if(b == sznan) b = a;
  if(a == sznan) panic("Unusable data.", __FILE__, __LINE__);

  /* Exact hit. */
  if(x[a] == x[b]) return x[a];

  /* Interpolate between two samples. */
  mdreal dB = (f[b] - q + 1e-10);
  mdreal dA = (q - f[a] + 1e-10);
  return (dB*(x[a]) + dA*(x[b]))/(dA + dB);
}
