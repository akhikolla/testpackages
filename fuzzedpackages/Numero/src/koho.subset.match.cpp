/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *  
 */
Point*
Subset::match(vector<Subset>& subsets, const vector<mdreal>& delta,
	      Point* pnt) {
  mdreal rlnan = medusa::rnan();
  
  /* Sort by distance. */
  vector<mdreal> tmp = delta;
  vector<mdsize> sorted = medusa::sortreal(tmp, 1);
  
  /* Find the best available subset. */
  for(vector<mdsize>::iterator it = sorted.begin();
      it != sorted.end(); it++) {
    if(*it >= subsets.size())
      panic("Invalid program state.", __FILE__, __LINE__);
    if(delta[*it] == rlnan) continue;
    Point* prev = subsets[*it].join(pnt, delta[*it]);
    if(prev != pnt) return prev;
  }
  return pnt;
}
