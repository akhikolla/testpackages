/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "punos.local.h"

/*
 *
 */
mdreal
Topology::sigma() const {
  TopologyBuffer* p = (TopologyBuffer*)buffer;
  return p->sigma;
}
