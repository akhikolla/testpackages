#include "utils_random_shuffle.h"

inline int RandWrapper(const int n) { return floor(unif_rand()*n); }
arma::vec RandomShuffle(arma::vec a) {
  arma::vec b = a;
  std::random_shuffle(b.begin(), b.end(), RandWrapper);
  return b;
}

