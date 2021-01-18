/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
mdreal
Trainer::cycle(vector<Point*>& points, const Topology& topo) {
  mdreal delta = this->match(points, topo);
  this->update(topo);
  return delta;
}
