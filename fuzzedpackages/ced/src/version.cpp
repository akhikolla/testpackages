#include <Rcpp.h>
#include "ced.h"

//' @title
//' Compact Encodig Detector Verion
//'
//' @description
//' Backed library version string.
//'
//' @return
//' Numeric version of the upstream library.
//'
//' @export
//'
// [[Rcpp::export(rng = false)]]
Rcpp::List ced_version() {
  const char* s = Version();
  std::stringstream ss(s);
  std::vector<int> tmp;
  std::string num;
  while (std::getline(ss, num, '.')) {
    tmp.push_back(std::stoi(num));
  }
  Rcpp::List res = Rcpp::List::create(tmp);
  res.attr("class") = "numeric_version";
  return res;
}
