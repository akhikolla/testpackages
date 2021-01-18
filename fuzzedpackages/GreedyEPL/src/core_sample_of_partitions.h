#ifndef SAMPLE_OF_PARTITIONS_H
#define SAMPLE_OF_PARTITIONS_H

// # include <armadillo>
# include <RcppArmadillo.h>

class sample_of_partitions {
public:
  sample_of_partitions() {};
  sample_of_partitions(arma::mat, arma::vec, arma::vec);

  unsigned int niter;
  unsigned int N;
  unsigned int Kup;
  arma::mat sample;
  arma::vec weights;
  double sum_of_weights;
  arma::mat sample_counts;
  arma::field<arma::vec> non_empty_groups_sample;
  arma::vec decision;
  arma::vec decision_counts;
  arma::vec non_empty_groups_decision;
  arma::cube contingency_tables;
  arma::vec losses;
  double epl_value;
  arma::vec deltas;
  std::string loss_function_name;

  void Print();
  void Summary();
  void EvaluateCounts();
  void EvaluateLosses();
  double EvaluateDelta(unsigned int, unsigned int);
  void EvaluateDeltas(unsigned int);
  void Move(unsigned int, unsigned int);
};


#endif // SAMPLE_OF_PARTITIONS_H
