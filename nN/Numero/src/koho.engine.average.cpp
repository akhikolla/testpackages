/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<vector<mdreal> >
Engine::average() const {
  EngineBuffer* p = (EngineBuffer*)buffer;
  mdreal rlnan = medusa::rnan();

  /* Check if data are ready. */
  if((p->cache).size() < 1) p->prepare();
  
  /* Estimate district averages for each column. */
  vector<mdsize>& loci = (p->bmus).first;
  vector<vector<mdreal> > planes((p->cache).size());
  for(mdsize j = 0; j < planes.size(); j++)
    planes[j] = (p->topology).diffuse(loci, p->cache[j]);

  /* Estimate smoothed point frequencies. */
  vector<vector<mdreal> > hgrams = this->histograms();
  if(planes.size() != hgrams.size())
    medusa::panic("Size mismatch.", __FILE__, __LINE__);
  
  /* Normalize values. */
  for(mdsize j = 0; j < planes.size(); j++) {
    vector<mdreal>& sums = planes[j];
    vector<mdreal>& nums = hgrams[j];
    for(mdsize i = 0; i < sums.size(); i++)
      if(nums[i] > 0.0) sums[i] /= (nums[i] + 1e-20);
      else sums[i] = rlnan;
  }
  return planes;
}
