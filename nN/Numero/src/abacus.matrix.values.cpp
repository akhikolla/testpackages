/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
void
Matrix::values(vector<Element>& elem) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  unordered_map<mdsize, Array>& rowdata = p->rowdata;
  for(mdsize k = 0; k < elem.size(); k++) {
    Element& e = elem[k];
    
     /* Check if symmetric. */
    mdsize r = e.row;
    mdsize c = e.column;
    if(p->symmflag && (r > c)) {
      r = e.column;
      c = e.row;
    }

    /* Check if outside bounds. */
    if(r >= p->nrows) {e.value = p->rlnan; continue;}
    if(c >= p->ncols) {e.value = p->rlnan; continue;}
    
    /* Check if row exists. */
    if(rowdata.count(r) < 1) {e.value = p->rlnan; continue;}
    
    /* Copy value. */
    e.value = rowdata[r][c];
  }
}
