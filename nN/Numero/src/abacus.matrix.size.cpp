/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdsize
Matrix::size() const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  return p->nrows;
}
