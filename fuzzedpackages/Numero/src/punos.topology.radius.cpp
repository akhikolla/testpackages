/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
mdreal
Topology::radius() const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  return p->maxradius;
}
