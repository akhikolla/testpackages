/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
Matrix::symmetric(const bool flag) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if((p->rowdata).size() < 1) p->symmflag = flag;
  else panic("Matrix not empty.", __FILE__, __LINE__);
}

/*
 *
 */
bool
Matrix::symmetric() const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  return p->symmflag;
}
