#include <Rcpp.h>
#include <string>
#include <vector>
#include "c2integer-utils.h"

using namespace Rcpp;

// Function to convert a Chinese numeral to integer
// [[Rcpp::export]]
long long c2integer_conv(const std::vector<std::string> number, const List conv_t)
{
  const DataFrame chr_t = as<DataFrame>(conv_t["chr_t"]);
  const DataFrame scale_t = as<DataFrame>(conv_t["scale_t"]);
  const std::string zero = conv_t["zero"];

  std::vector<long long> normal_stack {0}; // hold values before summing at the end
  std::vector<long long> scale_stack {0}; // hold values of premature summing for recursive scaling
  std::vector<int> scale_n; // the scale number of the char, 0 if not scale char
  int tmp = 0; // to hold value of number to be manipulated by scale digit
  std::string digit, scale_place;
  int previous_scale = 0;
  long long digit_n;

  for (unsigned int i = 0; i < number.size(); ++i) {
    digit = number[i];
    if (digit == zero) {
      // if numeral is zero
      scale_n.push_back(0);
      continue;
    }
    if (contains(chr_t["c"], digit)) {
      // if numeral is a digit
      digit = subset_df(chr_t, digit);
      scale_n.push_back(0);
      tmp = digit[0] - 48;
      continue;
    }
    // if numeral is a scale char
    scale_place = subset_df(scale_t, digit);
    scale_n.push_back(std::stoi(scale_place));
    digit_n = pow10(scale_n[i] - 1);
    if (i == 0) {
      // e.g.  十二
      normal_stack.push_back(digit_n);
      continue;
    }
    if (scale_n[i - 1] > 0) {
      // if previous numeral is also a scale char, e.g. 一百萬
      scale_stack.push_back(std::accumulate(normal_stack.begin(), normal_stack.end(), 0LL) * digit_n);
      normal_stack = {0};
    }
    if (!contains(scale_t["n"], scale_n[i] - 1)) {
      // if current scale is a large scale which requires recursive scaling, e.g. 一億 and 一兆
      scale_stack.push_back((std::accumulate(normal_stack.begin(), normal_stack.end(), 0LL) + tmp) * digit_n);
      normal_stack = {0};
      tmp = 0;
      continue;
    }
    if (i > 1) {
      for (auto j = scale_n.end() - 2; j >= scale_n.begin(); --j) { // find previous scale
        if (*j != 0) {
          previous_scale = *j;
          break;
        }
      }
      if (previous_scale < scale_n[i]) {
        // if previous scale is smaller than current scale, e.g. 一百二十三萬
        scale_stack.push_back((std::accumulate(normal_stack.begin(), normal_stack.end(), 0LL) + tmp) * digit_n);
        normal_stack = {0};
        tmp = 0;
        continue;
      }
    }
    // e.g. 一百二十三
    normal_stack.push_back(tmp * digit_n);
    tmp = 0;
  }
  return std::accumulate(normal_stack.begin(), normal_stack.end(), 0LL) +
    std::accumulate(scale_stack.begin(), scale_stack.end(), 0LL) + tmp;
}

// Function to convert a Chinese numeral to integer literally
// [[Rcpp::export]]
std::string c2integer_literal(const std::vector<std::string> number, const List conv_t)
{
  const DataFrame chr_t = as<DataFrame>(conv_t["chr_t"]);
  const std::string zero = conv_t["zero"];

  std::string converted;
  for (auto& n : number) {
    if (n == zero)
      converted += "0";
    else
      converted += subset_df(chr_t, n);
  }
  return converted;
}
