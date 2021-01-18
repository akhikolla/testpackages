/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
mdreal
Matrix::value(const mdsize r0, const mdsize c0) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;

  /* Check if symmetric. */
  mdsize r = r0;
  mdsize c = c0;
  if(p->symmflag && (r > c)) {
    r = c0;
    c = r0;
  }

  /* Check if outside bounds. */
  if(r >= p->nrows) return p->rlnan;
  if(c >= p->ncols) return p->rlnan;

  /* Check if row exists. */
  if((p->rowdata).count(r) < 1) return p->rlnan;

  /* Return value. */
  return p->rowdata[r][c];
}
