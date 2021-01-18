/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#include "abacus.local.h"

/*
 *
 */
class ValueComparator {
private:
  int flag;
public:
  ValueComparator() {this->flag = 0;};
  ValueComparator(const int k) {this->flag = k;};
  ~ValueComparator() {};
  bool operator()(const Element& x, const Element& y) {
    if(flag > 0) return (x.value < y.value);
    if(flag < 0) return (x.value > y.value);
    panic("Bad parameter.", __FILE__, __LINE__);
    return false;
  };
};

/*
 *
 */
vector<Element>
MatrixBuffer::elements(const int sortflag, const bool remflag) {
  
  /* Collect elements. */
  vector<Element> elem;
  for(unordered_map<mdsize, Array>::iterator it = rowdata.begin();
      it != rowdata.end(); it++) {
    (it->second).elements(elem, it->first);
    if(remflag) it->second = Array();
  }

  /* Clear contents. */
  if(remflag) {
    this->nrows = 0;
    this->ncols = 0;
    (this->rowdata).clear();
  }
  
  /* Sort according to value. */
  if(sortflag == 0) return elem;
  ValueComparator cmp(sortflag);
  sort(elem.begin(), elem.end(), cmp);
  return elem;
}
