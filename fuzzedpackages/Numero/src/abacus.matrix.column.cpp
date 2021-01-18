/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
vector<mdreal>
Matrix::column(const mdsize c) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  if(p->symmflag) panic("Symmetric matrix.", __FILE__, __LINE__);
  if(c >= p->ncols) return vector<mdreal>();
  vector<mdreal> output(p->nrows, p->rlnan);
  unordered_map<mdsize, Array>& rowdata = p->rowdata;
  for(unordered_map<mdsize, Array>::iterator it = rowdata.begin();
      it != rowdata.end(); it++)
    output[it->first] = it->second[c];
  return output;
}

/*
 *
 */
mdsize
Matrix::column(vector<Element>& output, const mdsize c) const {
  MatrixBuffer* p = (MatrixBuffer*)buffer;
  unordered_map<mdsize, Array>& rowdata = p->rowdata;
  mdreal rlnan = medusa::rnan();
  output.clear();
  if(p->symmflag) panic("Symmetric matrix.", __FILE__, __LINE__);
  if(c >= p->ncols) return 0;
  for(unordered_map<mdsize, Array>::iterator it = rowdata.begin();
      it != rowdata.end(); it++) {
    Element e;
    e.value = it->second[c];
    if(e.value == rlnan) continue;
    e.row = it->first;
    e.column = c;
    output.push_back(e);
  }
  return output.size();
}
