/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<mdreal>
Matrix::row(const mdsize r) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(p->symmflag) panic("Symmetric matrix.", __FILE__, __LINE__);
  if(r >= p->nrows) return vector<mdreal>();
  if((p->rowdata).count(r) < 1) return vector<mdreal>();
  vector<mdreal> output = (p->rowdata[r]).values();
  output.resize(p->ncols, p->rlnan);
  return output;
}

/*
 *
 */
mdsize
Matrix::row(vector<Element>& output, const mdsize r) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  output.clear();
  if(p->symmflag) panic("Symmetric matrix.", __FILE__, __LINE__);
  if(r >= p->nrows) return 0;
  if((p->rowdata).count(r) < 1) return 0;
  (p->rowdata[r]).elements(output, r);
  return output.size();
}
