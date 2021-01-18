/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
vector<mdsize>
Topology::neighbors(const mdsize ind) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  if(ind >= (p->network).size()) return vector<mdsize>();
  vector<mdsize> units;
  LinkMap::iterator pos;
  LinkMap neigh = p->network[ind];
  for(pos = neigh.begin(); pos != neigh.end(); pos++)
    units.push_back(pos->first);
  return units;
}
