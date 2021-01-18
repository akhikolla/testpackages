/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
void
Trainer::update(const Topology& topo) {
  mdreal rlnan = medusa::rnan();
  mdreal sigma = topo.sigma();
 
  /* Estimate centroids. */
  Matrix sums;
  Matrix numbers;
  Matrix centroids;
  for(mdsize i = 0; i < subsets.size(); i++) {
    vector<mdreal> mom1, mom2;
    vector<mdsize> nums = subsets[i].moments(mom1, mom2);
    for(mdsize j = 0; j < nums.size(); j++) {
      mdsize n = nums[j];
      if(n < 1) continue;
      mdreal m1 = mom1[j];
      mdreal mu = m1/n;
      sums.insert(i, j, m1);
      numbers.insert(i, j, 1.0*n);
      centroids.insert(i, j, mu);
    }
  }
  
  /* Apply spatial smoothing. */
  if(sigma > 0.0) {  
    vector<mdsize> bmus;
    for(mdsize i = 0; i < sums.size(); i++)
      bmus.push_back(i);
    
    /* Calculate smoothed component planes. */
    for(mdsize j = 0; j < sums.order(); j++) {
      vector<mdreal> sumplane = topo.diffuse(bmus, sums.column(j));
      vector<mdreal> numplane = topo.diffuse(bmus, numbers.column(j));
      
      /* Re-calculate mean values. */
      for(mdsize i = 0; i < sumplane.size(); i++) {
	if(numplane[i] == rlnan) continue;
	if(numplane[i] == 0.0) continue;
	mdreal value = (sumplane[i])/(numplane[i] + 1e-9);
	centroids.insert(i, j, value);
      }
    }
  }
  
  /* Update prototypes. */
  (this->prototypes).resize(subsets.size());
  for(mdsize i = 0; i < subsets.size(); i++) {
    vector<mdreal> mu = centroids.row(i);
    vector<mdreal>& proto = this->prototypes[i];
    if(proto.size() < mu.size()) proto.resize(mu.size(), rlnan);
    for(mdsize j = 0; j < mu.size(); j++)
      if(mu[j] != rlnan) proto[j] = mu[j];
  }
}
