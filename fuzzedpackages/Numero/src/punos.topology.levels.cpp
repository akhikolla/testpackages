/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
vector<mdreal>
Topology::levels() const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  return p->levels;
}
