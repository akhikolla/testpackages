/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
EmpiricalBuffer::contents(vector<mdreal>& x, vector<mdreal>& w) {

  /* Check if cached values are available. */
  if(valsorted.size() > 0) {
    x = valsorted;
    w = wsorted;
    return;
  }
  else {
    x.resize(data.size());
    w.resize(data.size());
  }
  
  /* Copy values. */
  mdsize ncopied = 0;
  unordered_map<mdreal, mdreal>::const_iterator pos;
  for(pos = data.begin(); pos != data.end(); pos++) {
    x[ncopied] = pos->first;
    w[ncopied] = pos->second;
    ncopied++;
  }

  /* Sort by value. */
  vector<mdreal> wcopy = w;
  vector<mdsize> sorted = sortreal(x, 1);
  for(mdreal k = 0; k < sorted.size(); k++)
    w[k] = wcopy[sorted[k]];

  /* Copy values to cache. */
  this->valsorted = x;
  this->wsorted = w;
}
