/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
mdreal
Topology::distance(const mdsize a, const mdsize b) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  vector<District>& coord = p->coord;
  mdsize ndistricts = coord.size();
  if(a >= ndistricts) return medusa::rnan();
  if(b >= ndistricts) return medusa::rnan();
  mdreal dx = (coord[b].x - coord[a].x);
  mdreal dy = (coord[b].y - coord[a].y);
  return sqrt(dx*dx + dy*dy);
}
