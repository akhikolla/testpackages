#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
CharacterVector ngrams(CharacterVector tokens, int n, std::string sep = "_")
{
  int len = tokens.size();
  CharacterVector out(len);
  int local_i = 0;

  for (int i = 0; i < len; i++)
  {
    if (i < n) {
      out[i] = CharacterVector::get_na();
    } else {
      out[i] = tokens[i];
      for (int j = 1; j < n; j++)
      {
        out[i] = std::string(tokens[i - j]) + sep + std::string(out[i]);
      }
    }
  }

  return out;
}
