/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
mdreal
Trainer::match(vector<Point*>& points, const Topology& topo) {
  mdreal rlnan = medusa::rnan();
  mdsize sznan = medusa::snan();

  /* Clear previous matches. */
  for(mdsize i = 0; i < subsets.size(); i++)
    subsets[i].clear();

  /* Reset cluster labels. */
  vector<Point*> pool;
  for(mdsize i = 0; i < points.size(); i++) {
    Point* pnt = points[i]; pnt->move(sznan);
    pool.push_back(pnt);
  }

  /* Check if prototypes are available. */
  if(prototypes.size() < 1) {
    mdsize nunits = topo.size();
    for(mdsize i = 0; i < pool.size(); i++)
      subsets[i%nunits].join(pool[i], 0.0);
    return rlnan;
  }
  
  /* Find the best available matches. */
  while(pool.size() > 0) {
    
    /* Process remaining points. */
    vector<Point*> rejected;
    for(mdsize i = 0; i < pool.size(); i++) {
      Point* pnt = pool[i];
      
      /* Estimate distances to centroids. */
      vector<mdreal> delta = this->distance(*pnt);

      /* Find the best available subset. */
      Point* replc = Subset::match(this->subsets, delta, pnt);
      if(replc == pnt) panic("Invalid program state.", __FILE__, __LINE__);
      if(replc != NULL) rejected.push_back(replc);
    }
    pool = rejected;
  }

  /* Estimate mean distance. */
  long double dsum = 0.0;
  long double wsum = 0.0;
  for(mdsize i = 0; i < points.size(); i++) {
    Point* pnt = points[i];
    mdsize bmu = pnt->location();
    mdreal d = this->distance(*pnt, bmu);
    if(d == rlnan) continue;
    dsum += d;
    wsum += 1.0;
  }
  if(wsum <= 0.0) return rlnan;
  return dsum/wsum;
}
