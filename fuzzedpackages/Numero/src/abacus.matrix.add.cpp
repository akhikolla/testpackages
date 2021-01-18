/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
bool
Matrix::add(const Element& e) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(e.value == p->rlnan) return false;

  /* Check if symmetric. */
  mdsize r = e.row;
  mdsize c = e.column;
  if(p->symmflag && (r > c)) {
    r = e.column;
    c = e.row;
  }
  
  /* Update contents. */
  (p->rowdata[r]).update(c, e.value, false);
  if(r >= p->nrows) p->nrows = (r + 1);
  if(c >= p->ncols) p->ncols = (c + 1);
  return true;
}



/*
 *
 */
bool
Matrix::add(const mdsize r0, const mdsize c0, const mdreal x) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(x == p->rlnan) return false;

  /* Check if symmetric. */
  mdsize r = r0;
  mdsize c = c0;
  if(p->symmflag && (r > c)) {
    r = c0;
    c = r0;
  }
  
  /* Update contents. */
  (p->rowdata[r]).update(c, x, false);
  if(r >= p->nrows) p->nrows = (r + 1);
  if(c >= p->ncols) p->ncols = (c + 1);
  return true;
}
