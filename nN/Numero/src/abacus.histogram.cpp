/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<mdreal>
abacus::histogram(const vector<mdreal>& x0, const vector<mdreal>& q0) {
  vector<mdreal> w0(x0.size(), 1.0);
  return histogram(x0, w0, q0);
}

/*
 *
 */
vector<mdreal>
abacus::histogram(const vector<mdreal>& x, const vector<mdreal>& w,
		  const vector<mdreal>& q) {
  mdsize nelem = x.size();
  mdsize nbins = q.size();
  mdsize sznan = medusa::snan();
  mdreal rlnan = medusa::rnan();
  vector<mdreal> h(nbins, 0.0);

  /* Check weights. */
  if(w.size() != nelem) {
    worry("Incompatible inputs.", __FILE__);
    return h;
  }

  /* Check if any data. */
  if(nelem < 1) return h;  
  if(nbins < 2) return h;

  /* Check that bins are sorted and distinct. */
  for(mdsize k = 1; k < nbins; k++) {
    if(q[k] > q[k-1]) continue;
    worry("Unusable bin positions.", __FILE__);
    return h;
  }

  /* Collect weighted hit counts. */
  vector<mdreal> hits(nbins, 0.0);
  for(mdsize i = 0; i < nelem; i++) {
    if(x[i] == rlnan) continue;
    long double key = x[i];

    /* Find bin edges. */
    Site slot = binsearch(q, key);

    /* Collapse "tails" to nearest bin edge. */
    mdsize a = slot.bounds.first;
    mdsize b = slot.bounds.second;
    if(a == sznan) a = b;
    if(b == sznan) b = a;
    if(a == sznan) continue;

    /* Exact hit. */
    if(a == b) {
      hits[a] += w[i];
      continue;
    }
    
    /* Interpolate hit contribution. */
    double dA = (key - q[a] + 1e-10);
    double dB = (q[b] - key + 1e-10);
    hits[a] += dB*(w[i])/(dA + dB);
    hits[b] += dA*(w[i])/(dA + dB);
  }
  return hits;
}
