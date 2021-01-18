/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<vector<mdreal> >
Engine::histograms() const {
  EngineBuffer* p = (EngineBuffer*)buffer;
  mdreal rlnan = medusa::rnan();

  /* Check if data are ready. */
  if((p->cache).size() < 1) p->prepare();
  if((p->freqs).size() > 0) return p->freqs;

  /* Estimate point frequences for each data column. */
  vector<mdsize>& loci = (p->bmus).first;
  for(mdsize j = 0; j < (p->cache).size(); j++) {
    vector<mdreal> ones = p->cache[j]; /* deep copy */

    /* Indicators for usable values. */
    for(mdsize k = 0; k < ones.size(); k++)
      ones[k] = (ones[k] != rlnan);

    /* Smoothed point histogram on the map. */
    vector<mdreal> h = (p->topology).diffuse(loci, ones);

    /* If dataset is complete, one histogram is enough. */
    if(p->complete) {
      p->freqs.resize((p->cache).size(), h);
      return p->freqs;
    }

    /* Add new histogram. */
    (p->freqs).push_back(h);
  }
  return p->freqs;
}
