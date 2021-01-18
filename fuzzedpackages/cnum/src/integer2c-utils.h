#include <Rcpp.h>
#include <string>
#include <boost/locale/encoding_utf.hpp>
#include "utils.h"

// [[Rcpp::depends(BH)]]

// Function to convert string to wstring
std::wstring s2ws(const std::string& str)
{
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<wchar_t>(str.c_str(), str.c_str() + str.size());
}

// Function to convert wstring to string
std::string ws2s(const std::wstring& wstr)
{
  using boost::locale::conv::utf_to_utf;
  return utf_to_utf<char>(wstr.c_str(), wstr.c_str() + wstr.size());
}

// Function to return index of a matching element in a numeric vector
int subset_num(const Rcpp::NumericVector x, const int y)
{
  int index = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == y + 1) {
      index = i;
      break;
    }
  }
  return index;
}

// Function to return the corresponding chr of df$c given n of df$n
std::string subset_df(const Rcpp::DataFrame df, const int n)
{
  const Rcpp::NumericVector n_col = df["n"];
  const Rcpp::CharacterVector c_col = df["c"];
  int index = subset_num(n_col, n);
  return Rcpp::as<std::string>(c_col[index]);
}
