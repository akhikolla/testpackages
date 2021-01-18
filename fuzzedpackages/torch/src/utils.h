#include "torch_types.h"

template <class type>
Rcpp::XPtr<type> make_xptr  (type x) {
  auto * out = new type(x);
  return Rcpp::XPtr<type>(out);
}

template <class type>
Rcpp::XPtr<type> make_xptr  (type x, std::string dyn_type) {
  auto * out = new type(x);
  auto ptr = Rcpp::XPtr<type>(out);
  ptr.attr("dynamic_type") = dyn_type;
  return ptr;
}

template <class type, int n>
std::array<type, n> std_vector_to_std_array (std::vector<type> x) {
  std::array<type,n> out;
  std::copy_n(x.begin(), n, out.begin());
  return out;
}
