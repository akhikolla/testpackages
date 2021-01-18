#include <Rcpp.h>
#include <string>
#include "utils.h"

// Function to calulate 10 to the power of given exponent
long long pow10(const int exp)
{
  long long result = 1;
  for (int i = 1; i <= exp; ++i)
    result *= 10;
  return result;
}

// Function to return index of a matching element in a character vector
int subset_chr(const Rcpp::CharacterVector x, const std::string y)
{
  int index = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] == y) {
      index = i;
      break;
    }
  }
  return index;
}

// Function to return the corresponding chr of df$n given c of df$c
std::string subset_df(const Rcpp::DataFrame df, const std::string c)
{
  const Rcpp::NumericVector n_col = df["n"];
  const Rcpp::CharacterVector c_col = df["c"];
  if (!contains(c_col, c))
    return "";
  int index = subset_chr(c_col, c);
  int output = n_col[index];
  return std::to_string(output);
}
