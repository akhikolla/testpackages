/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
mdsize
Topology::depth() const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  return (p->levels).size();
}
