#ifndef NORMALISED_VARIATION_OF_INFORMATION_H
#define NORMALISED_VARIATION_OF_INFORMATION_H

// # include <armadillo>
# include <RcppArmadillo.h>
# include "core_sample_of_partitions.h"
# include "utils_entropy.h"

class normalised_variation_of_information : public sample_of_partitions {
public:
  normalised_variation_of_information() {};
  normalised_variation_of_information(arma::mat, arma::vec, arma::vec);

  double entropy_decision;
  arma::vec entropies_sample;
  arma::vec joint_entropies;

  void EvaluateLosses();
  double EvaluateDelta(unsigned int, unsigned int);
  void EvaluateDeltas(unsigned int);
  void Move(unsigned int, unsigned int);
};

#endif // NORMALISED_VARIATION_OF_INFORMATION_H
