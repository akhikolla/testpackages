/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
vector<mdreal>
Model::prototype(const mdsize unit) const {
  ModelBuffer* p = (ModelBuffer*)buffer;
  return (p->codebook).row(unit);
}
