/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
Site
Topology::stratify(const mdreal z) const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  return medusa::binsearch(p->levels, z);
}
