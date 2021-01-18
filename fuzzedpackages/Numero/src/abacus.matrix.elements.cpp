/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<Element>
Matrix::elements(const int flag) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  return p->elements(flag, false);
}
