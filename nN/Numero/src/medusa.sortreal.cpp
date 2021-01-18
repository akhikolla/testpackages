/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
struct RealItem {
  bool isnan;
  mdreal value;
  mdsize rank;
};

/*
 *
 */
class RealCompare {
public:
  bool flag;
public:
  bool operator()(const RealItem& x, const RealItem& y) const {
    if(x.isnan != y.isnan) return y.isnan;
    if(flag) return (x.value < y.value);
    return (x.value > y.value);
  };
};

/*
 *
 */
vector<mdsize>
medusa::sortreal(vector<mdreal>& items, const int direction) {
  if(direction == 0) panic("Unusable input.", __FILE__, __LINE__);
  mdsize nelem = items.size();
  mdreal rlnan = medusa::rnan();

  /* Count sorted items. */
  mdsize nsorted = 1;
  RealCompare cmp; cmp.flag = (direction > 0);
  for(mdsize k = 1; k < nelem; k++) {
    mdreal xA = items[k-1];
    mdreal xB = items[k];
    if(cmp.flag) {if(xA > xB) break;}
    else {if(xA < xB) break;}
    nsorted++;
  }

  /* Check if already sorted. */
  vector<mdsize> array; array.reserve(nelem); 
  if(nsorted >= nelem) {
    for(mdsize k = 0; k < nelem; k++)
      array.push_back(k);
    return array;
  }

  /* Create float-rank pairs. */
  vector<RealItem> tuples(nelem);
  for(mdsize k = 0; k < nelem; k++) {
    tuples[k].isnan = (items[k] == rlnan);
    tuples[k].value = items[k];
    tuples[k].rank = k;
  }

  /* Sort items. */
  sort(tuples.begin(), tuples.end(), cmp);

  /* Collect indices. */
  for(mdsize k = 0; k < nelem; k++) {
    array.push_back(tuples[k].rank);
    items[k] = tuples[k].value;
  }
  return array;
}
