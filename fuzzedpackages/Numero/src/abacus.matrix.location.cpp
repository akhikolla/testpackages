/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdsize
Matrix::location(const string& key, const string& flag) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(flag == "row") return (p->rownames).rank(key);
  if(flag == "column") return (p->colnames).rank(key);
  return medusa::snan();
}
