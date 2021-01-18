#include <Rcpp.h>
using namespace Rcpp;

//' Substrings of a Character Vector (C++ Style)
//' 
//' @description Extract substrings from a \code{character} vector in C++.
//' 
//' @param x A \code{character} vector.
//' @param pos The start point of the substring as \code{integer}. Position 
//' indications start from \code{1}, which is the default in R.
//' @param len The length of the substring as \code{integer}.
//' 
//' @return
//' A \code{character} vector of the same length as 'x'.
//' 
//' @seealso
//' \url{http://www.cplusplus.com/reference/string/string/substr/},
//' \code{\link{substr}}.
//' 
//' @examples
//' substrC("Hello, world!", pos = 1, len = 5)
//' 
//' @export
// [[Rcpp::export]]
Rcpp::StringVector substrC(Rcpp::StringVector x, int pos, int len) {

  // translate position indication from R to C++ 
  pos = pos - 1;
  
  // initialize variables
  std::string s;

  int n = x.size();
  Rcpp::StringVector out(n);
  
  // cycle through elements in 'x'
  for (int i = 0; i < n; i++) {
    s = as<std::string>(x[i]);
    out[i] = s.substr(pos, len);
  }

  return out;
}
