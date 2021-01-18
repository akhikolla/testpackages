/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

#define ATOMSIZE 8

static mdsize binsearch_exec(const vector<mdreal>*, const mdreal,
			     const mdsize, const mdsize);


/*
 * Input must be sorted in ascending order, it will not be checked
 * due to performance penalty!
 */
Site
medusa::binsearch(const vector<mdreal>& x, const mdreal key) {
  mdsize nelem = x.size();
  mdsize sznan = medusa::snan();
  mdreal rlnan = medusa::rnan();

  /* Empty locus. */
  Site slot;
  slot.bounds.first = sznan;
  slot.bounds.second = sznan;
  slot.weights.first = rlnan;
  slot.weights.second = rlnan;
  slot.usable = 0;

  /* Check if search is needed. */
  if(nelem < 1) return slot;
  if(key == rlnan) return slot;

  /* Find adjacent bin edge. */
  mdsize pos = binsearch_exec(&x, key, 0, (nelem - 1));

  /* Exact hit. */
  if(key == x[pos]) {
    slot.bounds.first = pos;
    slot.bounds.second = pos;
    slot.weights.first = 1.0;
    slot.weights.second = 0.0;
    slot.usable = 2;
    return slot;
  }

  /* Smaller edge. */
  mdsize& a = slot.bounds.first;
  mdsize& b = slot.bounds.second;
  if(key > x[pos]) {
    a = pos;
    b = (pos + 1);
    if(b >= nelem) b = sznan;
  }

  /* Larger edge. */
  if(key < x[pos]) {
    if(pos < 1) a = sznan;
    else a = (pos - 1);
    b = pos;
  }

  /* Check if duplicates below. */
  if(a != sznan) {
    while(x[a] == x[pos]) {
      if(a < 1) break;
      if(x[a-1] != x[pos]) break;
      a--;
    }
    if(x[a] == rlnan) a = sznan;
    else slot.usable += 1;
  }

  /* Check if duplicates above. */
  if(b != sznan) {
    while(x[b] == x[pos]) {
      if((b + 1) >= nelem) break;
      if(x[b+1] != x[pos]) break;
      b++;
    }
    if(x[b] == rlnan) b = sznan;
    else slot.usable += 1;
  }

  /* Check if enough usable values. */
  if(slot.usable < 2) return slot;

  /* Estimate interpolation weights. */
  double dA = (key - x[a]);
  double dB = (x[b] - key);
  slot.weights.first = dB/(dA + dB);
  slot.weights.second = dA/(dA + dB);
  return slot;
}

/*
 *
 */
mdsize
binsearch_exec(const vector<mdreal>* p_q, const mdreal key,
	       const mdsize alpha, const mdsize omega) {

  /* Check if end of recursion. */  
  if((omega - alpha) < ATOMSIZE) {
    for(mdsize k = alpha; k <= omega; k++) {
      if(key > (*p_q)[k]) continue;
      return k;
    }
    return omega;
  }

  /* Determine local center. */
  mdsize pivot = (alpha + omega)/2;

  /* Split search region. */
  mdreal x = (*p_q)[pivot];
  if(key < x) return binsearch_exec(p_q, key, alpha, pivot);
  if(key > x) return binsearch_exec(p_q, key, pivot, omega);
  return pivot;
}
