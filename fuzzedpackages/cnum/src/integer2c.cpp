#include <Rcpp.h>
#include <string>
#include <regex>
#include "integer2c-utils.h"

using namespace Rcpp;

// Function to convert integer to Chinese numeral
// [[Rcpp::export]]
std::string integer2c(const std::string number, const List conv_t)
{
  const std::wstring zero = s2ws(conv_t["zero"]);

  if (number == "0")
    return ws2s(zero);

  const DataFrame chr_t = as<DataFrame>(conv_t["chr_t"]);
  const DataFrame scale_t = as<DataFrame>(conv_t["scale_t"]);
  const std::wstring one = s2ws(subset_df(chr_t, 0));
  const std::wstring ten = s2ws(subset_df(scale_t, 1));

  const int n = number.size() - 1;
  std::wstring converted, sub_converted;
  int i = 0;
  int scale_n = n;
  int digit, lower_scale, scale_diff, sub_number_n;
  NumericVector n_col = scale_t["n"], lowers;
  std::string sub_number;

  while (i <= n) {
    if (contains(scale_t["n"], scale_n + 1)) {
      // non-recursive
      digit = number[i] - 48;
      if (digit == 0)
        converted += zero;
      else {
        converted += s2ws(subset_df(chr_t, digit - 1)); // respective Chinese numeral of the digit
        converted += s2ws(subset_df(scale_t, scale_n)); // respective Chinese numeral of the scale
      }
      ++i;
      --scale_n;
    } else {
      // recursive, e.g. process the bracketed digits (100)0000 to form 百萬
      lowers = n_col[(scale_n + 1) > n_col]; // find all lower scales
      lower_scale = max(lowers) - 1; // the next lower scale
      scale_diff = scale_n - lower_scale; // difference of current and next lower scale
      sub_number = number.substr(i, scale_diff + 1); // sub-number for recursive conversion
      sub_number_n = std::stoi(sub_number);
      if (sub_number_n > 0) {
        sub_converted = s2ws(integer2c(sub_number, conv_t)); // recursive_conversion
        if ((10 <= sub_number_n) && (sub_number_n <= 19))
          // add one when neccessary, e.g. 一兆零(一)十三億
          sub_converted = one + sub_converted;
        if (std::regex_search(sub_number, std::regex("^0")))
          // add zero for scales of digit zero, e.g. 一兆(零)三十三億
          sub_converted = zero + sub_converted;
        converted = (converted + sub_converted) + s2ws(subset_df(scale_t, lower_scale));
      }
      i += scale_diff + 1;
      scale_n -= scale_diff + 1;
    }
  }
  // correct "一十一" to "十一"
  converted = std::regex_replace(converted, std::wregex((L"^" + one) + ten), ten);
  // remove redundant zeros
  converted = std::regex_replace(converted, std::wregex(zero + L"{2,}"), zero);
  converted = std::regex_replace(converted, std::wregex(zero + L"$"), L"");
  return ws2s(converted);
}

// Function to convert integer to Chinese numeral literally
// [[Rcpp::export]]
std::string integer2c_literal(const std::string number, const List conv_t)
{
  const DataFrame chr_t = as<DataFrame>(conv_t["chr_t"]);
  const std::string zero = conv_t["zero"];

  std::string converted;
  int digit;
  for (auto& n : number) {
    digit = n - 48;
    if (digit == 0)
      converted += zero;
    else
      converted += subset_df(chr_t, digit - 1); // respective Chinese numeral of the digit
  }
  return converted;
}
