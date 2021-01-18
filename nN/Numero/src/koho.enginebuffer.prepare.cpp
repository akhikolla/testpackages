/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
void
EngineBuffer::prepare() {
  
  /* Discard previous contents. */
  (this->bmus).first.clear();
  (this->bmus).second.clear();
  (this->freqs).clear();
  (this->cache).clear();

  /* Allocate column vectors. */
  (this->cache).resize(order);
   
  /* Collect data. */
  vector<mdsize> loci;
  map<string, Point>::iterator pt;
  for(pt = points.begin(); pt != points.end(); pt++) {
    
    /* Copy location. */
    mdsize bmu = (pt->second).location();
    if(bmu >= topology.size()) continue; 
    loci.push_back(bmu);

    /* Copy values. */
    vector<mdreal> array = (pt->second).data();
    if(array.size() != order)
      medusa::panic("Unusable data point.", __FILE__, __LINE__);
    for(mdsize j = 0; j < array.size(); j++)
      (this->cache[j]).push_back(array[j]);
  }

  /* Prepare shuffling mask. */
  (this->bmus).first = loci;
  (this->bmus).second = loci;
}
