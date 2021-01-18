/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
Matrix::rename(const mdsize pos, const string& name, const string& flag) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(flag == "row") {
    if(name.size() < 1) (p->rownames).erase(pos);
    else (p->rownames).insert(pos, name);
  }
  if(flag == "column") {
    if(name.size() < 1) (p->colnames).erase(pos);
    else (p->colnames).insert(pos, name);
  }
}
