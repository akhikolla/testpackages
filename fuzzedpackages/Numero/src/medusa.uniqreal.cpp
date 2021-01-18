/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "medusa.local.h"

/*
 *
 */
vector<mdsize>
medusa::uniqreal(vector<mdreal>& items) {
  mdreal rlnan = medusa::rnan();
  vector<mdsize> mask;
  unordered_set<mdreal> bitset;
  for(mdsize i = 0; i < items.size(); i++) {
    if(items[i] == rlnan) continue;
    if(bitset.count(items[i]) > 0) continue;
    items[mask.size()] = items[i];
    bitset.insert(items[i]);
    mask.push_back(i);
  }
  items.resize(mask.size());
  return mask;
}
