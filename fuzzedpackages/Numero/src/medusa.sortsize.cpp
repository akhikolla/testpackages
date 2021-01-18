/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
struct SizeItem {
  bool isnan;
  mdsize value;
  mdsize rank;
};

/*
 *
 */
class SizeCompare {
public:
  bool flag;
public:
  bool operator()(const SizeItem& x, const SizeItem& y) const {
    if(x.isnan != y.isnan) return y.isnan;
    if(flag) return (x.value < y.value);
    return (x.value > y.value);
  };
};

/*
 *
 */
vector<mdsize>
medusa::sortsize(vector<mdsize>& items, const int direction) {
  if(direction == 0) panic("Unusable input.", __FILE__, __LINE__);  
  mdsize sznan = medusa::snan();

  /* Create integer-rank pairs. */ 
  mdsize nelem = items.size();
  vector<SizeItem> tuples(nelem);
  for(mdsize k = 0; k < nelem; k++) {
    tuples[k].isnan = (items[k] == sznan);
    tuples[k].value = items[k];
    tuples[k].rank = k;
  }
  
  /* Sort items. */
  SizeCompare cmp; cmp.flag = (direction > 0); 
  sort(tuples.begin(), tuples.end(), cmp);

  /* Collect indices. */
  vector<mdsize> mask(nelem); 
  for(mdsize k = 0; k < nelem; k++) {
    mask[k] = tuples[k].rank;
    items[k] = tuples[k].value;
  }
  return mask;
}
