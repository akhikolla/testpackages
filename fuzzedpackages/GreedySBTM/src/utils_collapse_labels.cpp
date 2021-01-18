#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List cpp_CollapseLabels (arma::vec vec) {
  arma::vec decision = vec;
  unsigned int i;
  unsigned int N = decision.size();
  unsigned int K = max(decision) + 1;
  for (i=0; i<N; ++i) if (decision.at(i) < 0) throw std::runtime_error("Decision vector provided has negative entries");
  arma::vec map_vector;
  map_vector.set_size(K);
  map_vector.fill(K);
  arma::vec decision_new = decision;
  unsigned int current_group_index = 0;
  for (i=0; i<N; ++i)
  {
    if (map_vector.at(decision.at(i)) < K) decision_new.at(i) = map_vector.at(decision.at(i));
    else
    {
      map_vector.at(decision.at(i)) = current_group_index;
      decision_new.at(i) = current_group_index;
      current_group_index += 1;
    }
  }
  return (Rcpp::List::create(Rcpp::Named("vec") = decision_new));
}

