#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

int longest_list_item(List l){
  int max = 0;
  for (int i = 0; i < l.size(); i++) {
    CharacterVector v = as<CharacterVector>(l[i]);
    if (max < v.size()) max = v.size();
  }
  return max;
}

// [[Rcpp::export]]
CharacterVector no_na_paste(List l, std::string sep = " "){
  int max = longest_list_item(l);
  CharacterVector out(max);
  CharacterVector v;
  std::string add;

  for (int li = 0; li < l.size(); li++) {
    v = as<CharacterVector>(l[li]);
    if (v.size() == 1) v = CharacterVector(max, v[0]);
    if (v.size() != max) ::Rf_error("List elements should have same length or length 1 (which is repeated)");
    for (int vi = 0; vi < max; vi++) {
      add = std::string(v[vi]);
      if (add != "NA") {
        if (out[vi] != "") add = sep + add;
        out[vi] = std::string(out[vi]) + add;
      }
    }
  }
  return out;
}
