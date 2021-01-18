#ifndef VARIATION_OF_INFORMATION_H
#define VARIATION_OF_INFORMATION_H

// # include <armadillo>
# include <RcppArmadillo.h>
# include "utils_entropy.h"
# include "core_sample_of_partitions.h"

class variation_of_information : public sample_of_partitions {
public:
  variation_of_information() {};
  variation_of_information(arma::mat, arma::vec, arma::vec);
  void EvaluateLosses();
  double EvaluateDelta(unsigned int, unsigned int);
  void EvaluateDeltas(unsigned int);
  void Move(unsigned int, unsigned int);
};

#endif // VARIATION_OF_INFORMATION_H
