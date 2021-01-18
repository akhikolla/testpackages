#include <Rcpp.h>

// Function to return whether y is an element of x
template <typename T>
bool contains(const Rcpp::CharacterVector x, const T y)
{
  return std::find(x.begin(), x.end(), y) != x.end();
}
