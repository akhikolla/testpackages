/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<mdsize>
Subset::moments(vector<mdreal>& sums1, vector<mdreal>& sums2) const {
  mdreal rlnan = medusa::rnan();
  sums1.clear();
  sums2.clear();
  
  /* Calculate base statistics. */
  vector<mdsize> nums;
  for(ContentMap::const_iterator it = contents.begin();
      it != contents.end(); it++) {
    const vector<Point*>& batch = it->second;
    for(mdsize i = 0; i < batch.size(); i++) {
      vector<mdreal> values = batch[i]->data();
      
      /* Check capacity. */
      mdsize ndim = values.size();
      if(nums.size() < ndim) {
	sums1.resize(ndim, 0.0);
	sums2.resize(ndim, 0.0);
	nums.resize(ndim, 0);
      }
      
      /* Update statistics. */
      for(mdsize j = 0; j < ndim; j++) {
	mdreal x = values[j];
	if(x == rlnan) continue;
	sums1[j] += x;
	sums2[j] += x*x;
	nums[j] += 1;
      }
    }
  }
  return nums;
}
