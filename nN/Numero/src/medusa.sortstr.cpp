/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
struct StringItem {
  double number;
  string value;
  mdsize rank;
};

/*
 *
 */
class StringCompare {
public:
  bool flag;
public:
  bool operator()(const StringItem& x, const StringItem& y) const {
    if(x.number != y.number) {
      if(flag) return (x.number < y.number);
      return (x.number > y.number);
    }
    if(flag) return (x.value < y.value);
    return (x.value > y.value);
  };
};

/*
 *
 */
vector<mdsize>
medusa::sortstr(vector<string>& items, const int direction) {
  if(direction == 0) panic("Unusable input.", __FILE__, __LINE__);  
  
  /* Create integer-rank pairs. */
  mdsize nelem = items.size();
  vector<StringItem> tuples(nelem);
  for(mdsize k = 0; k < nelem; k++) {
    tuples[k].number = atof(items[k].c_str());
    tuples[k].value = items[k];
    tuples[k].rank = k;
  }
  
  /* Sort items. */
  StringCompare cmp; cmp.flag = (direction > 0); 
  sort(tuples.begin(), tuples.end(), cmp);

  /* Collect indices. */
  vector<mdsize> mask(nelem); 
  for(mdsize k = 0; k < nelem; k++) {
    mask[k] = tuples[k].rank;
    items[k] = tuples[k].value;
  }
  return mask;
}
