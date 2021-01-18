#include <Rcpp.h>
using namespace Rcpp;

// mnis_bom
//' Strip out BOM from JSON data
//'
//' @param x The GET return to strip BOM out of
//' @export
// [[Rcpp::export]]
std::string mnis_bom(std::string x) {
  if (x.size() < 3)
    return x;

  if (x[0] == '\xEF' && x[1] == '\xBB' && x[2] == '\xBF')
    return x.substr(3);

  return x;
}

/*** R
x <- "\uFEFFabcdef"
print(x)
  print(mnis_bom(x))
  identical(x, mnis_bom(x))
  utf8ToInt(x)
  utf8ToInt(mnis_bom(x))
  */
