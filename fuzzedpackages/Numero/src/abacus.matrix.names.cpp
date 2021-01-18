/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<string>
Matrix::names(const string& flag) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  vector<string> output;
  if(flag == "row") {
    TwowayMap& rownames = p->rownames;
    for(mdsize i = 0; i < p->nrows; i++)
      output.push_back(rownames.name(i));
  }
  if(flag == "column") {
    TwowayMap& colnames = p->colnames;
    for(mdsize j = 0; j < p->ncols; j++)
      output.push_back(colnames.name(j));
  }
  return output;
}
