#include <Rcpp.h>
#include <random>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// returns the specified quantile, taking into account the weights of cusum values
inline double quantile_impl(double q /* the specified quantile, e.g. 0.5 = median */,
                     double total_weight /* the sum of the weights in sorted_c */,
                     std::vector<std::pair<double, double> > sorted_c /* sorted list of cusum values with weights */) {
  double target_value = q * total_weight;
  double accum = 0.;
  std::size_t i = 0;
  while(accum + sorted_c[i].second < target_value) {
    accum += sorted_c[i].second;
    i++;
  }
  return sorted_c[i].first;
}

// Given a list of cusum values with weights and a vector of
// required quantiles, returns a vector of those quantiles
inline NumericVector quantile(std::vector<std::pair<double, double> >& unsorted_c,
                       NumericVector& qs) {
  std::sort(unsorted_c.begin(), unsorted_c.end(),
            [](std::pair<double, double> const& a, std::pair<double, double> const& b) { return a.first < b.first; });
  NumericVector res(qs.size());
  double total_weight = 0.;
  for(auto const& elem : unsorted_c) { total_weight += elem.second; }
  for(int i = 0; i != qs.size(); ++i) {
    res[i] = quantile_impl(qs[i], total_weight, unsorted_c);
  }
  return res;
}
