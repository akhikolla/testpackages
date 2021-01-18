/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
vector<mdsize>
medusa::uniqstr(vector<string>& items) {
  vector<mdsize> mask;
  unordered_set<string> bitset;
  for(mdsize i = 0; i < items.size(); i++) {
    if(bitset.count(items[i]) > 0) continue;
    items[mask.size()] = items[i];
    bitset.insert(items[i]);
    mask.push_back(i);
  }
  items.resize(mask.size());
  return mask;
}
