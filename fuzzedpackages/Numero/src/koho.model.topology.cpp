/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Topology
Model::topology() const {
  ModelBuffer* p = (ModelBuffer*)buffer;
  return p->topology;
}
