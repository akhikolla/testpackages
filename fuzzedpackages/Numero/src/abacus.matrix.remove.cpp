/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<Element>
Matrix::remove(const int flag) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  return p->elements(flag, true);
}
  
/*
 *
 */
mdreal
Matrix::remove(const mdsize r0, const mdsize c0) {
  MatrixBuffer* p = (MatrixBuffer*)buffer;

  /* Check if symmetric. */
  mdsize r = r0;
  mdsize c = c0;
  if(p->symmflag && (r > c)) {
    r = c0;
    c = r0;
  }
  
  /* Check if row exists. */
  unordered_map<mdsize, Array>& rowdata = p->rowdata;
  if(rowdata.count(r) < 1) return p->rlnan;

  /* Check if column exists. */
  Array& array = p->rowdata[r];
  if(array.length() <= c) return p->rlnan; 

  /* Erase old value. */
  mdreal value = array.remove(c);
  
  /* Re-calculate number of rows. */
  if(array.length() < 1) {
    rowdata.erase(r); p->nrows = 0;
    for(unordered_map<mdsize, Array>::iterator it = rowdata.begin();
	it != rowdata.end(); it++)
      if(it->first >= p->nrows) p->nrows = (it->first + 1);
  }

  /* Re-calculate number of columns. */
  if(c == (p->ncols - 1)) {
    p->ncols = 0;
    for(unordered_map<mdsize, Array>::iterator it = rowdata.begin();
	it != rowdata.end(); it++) {
      mdsize n = (it->second).length();
      if(n > p->ncols) p->ncols = n;
    }
  }
  return value;
}
