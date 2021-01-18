#ifndef BINDER_H
#define BINDER_H

// # include <armadillo>
# include <RcppArmadillo.h>
# include "core_sample_of_partitions.h"

class binder : public sample_of_partitions {
public:
  binder() {};
  binder(arma::mat, arma::vec, arma::vec);

  void EvaluateLosses();
  double EvaluateDelta(unsigned int, unsigned int);
  void EvaluateDeltas(unsigned int);
  void Move(unsigned int, unsigned int);
};


#endif // BINDER_H
