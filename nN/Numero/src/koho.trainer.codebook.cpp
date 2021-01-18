/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "koho.local.h"

/*
 *
 */
Matrix
Trainer::codebook() const {
  Matrix mtx;
  for(mdsize i = 0; i < prototypes.size(); i++) {
    const vector<mdreal>& values = prototypes[i];
    for(mdsize j = 0; j < values.size(); j++)
      mtx.insert(i, j, values[j]);
  }
  return mtx;
}
