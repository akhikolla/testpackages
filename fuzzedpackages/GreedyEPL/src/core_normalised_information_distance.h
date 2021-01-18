#ifndef NORMALISED_INFORMATION_DISTANCE_H
#define NORMALISED_INFORMATION_DISTANCE_H

// # include <armadillo>
# include <RcppArmadillo.h>
# include "core_sample_of_partitions.h"
# include "utils_entropy.h"

class normalised_information_distance : public sample_of_partitions {
public:
  normalised_information_distance() {};
  normalised_information_distance(arma::mat, arma::vec, arma::vec);

  double entropy_decision;
  arma::vec entropies_sample;
  arma::vec joint_entropies;

  void EvaluateLosses();
  double EvaluateDelta(unsigned int, unsigned int);
  void EvaluateDeltas(unsigned int);
  void Move(unsigned int, unsigned int);
};

#endif // NORMALISED_INFORMATION_DISTANCE_H
