/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
mdreal
Topology::weight(const mdsize a, const mdsize b) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  mdsize ncoord = (p->coord).size();
  if(a >= ncoord) return medusa::rnan();
  if(a >= (p->network).size()) return medusa::rnan();
  LinkMap& neigh = p->network[a];
  LinkMap::iterator pos = neigh.find(b);
  if(pos == neigh.end()) return 0.0;
  return pos->second;
}
