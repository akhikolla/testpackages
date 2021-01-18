/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
bool
Engine::shuffle(const bool flag) {
  EngineBuffer* p = (EngineBuffer*)buffer;
  mt19937& twister = p->twister;

  /* Check engine state. */
  vector<mdsize>& loci = (p->bmus).first;
  mdsize nelem = loci.size();
  if(nelem < 1) return false;
  
  /* If some data are missing or if sampling with replacement,
     histograms need to be updated. */
  if(p->complete == false) (p->freqs).clear();
  if(flag == true) (p->freqs).clear();
  
  /* Shuffle point locations without replacement. */
  if(flag == false) {
    for(mdsize i = 0; i < nelem; i++) {
      mdsize ind = twister()%nelem;
      mdsize rank = loci[ind];
      loci[ind] = loci[i];
      loci[i] = rank;
    }
    return true;
  }

  /* Sampling with replacement. */
  const vector<mdsize>& orig = (p->bmus).second;
  for(mdsize i = 0; i < nelem; i++)
    loci[i] = orig[twister()%nelem];
  return true;
}
