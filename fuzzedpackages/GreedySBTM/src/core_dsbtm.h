#ifndef DSBTM_H
#define DSBTM_H

#include <RcppArmadillo.h>
#include "utils_bernoulli_marginal_delta.h"
#include "utils_random_shuffle.h"

class dsbtm
{
public:
  dsbtm();
  dsbtm(arma::cube, arma::mat, unsigned int, bool);
  
  // Data
  arma::cube adj;
  arma::cube active;
  unsigned int T;
  unsigned int N;
  
  // Allocations
  unsigned int Kup;
  arma::mat z;
  arma::mat counts;
  arma::vec total_counts;
  arma::vec non_empty_groups, non_empty_groups_no_0;
  arma::mat transition_counts;
  arma::vec transition_total_counts;
  unsigned int K, K_no_0;
  
  // Likelihood
  arma::mat eta;// SBM part: counts the number of edges between any two groups
  arma::mat zeta;// SBM part: counts the number of non-edges between any two groups
  arma::mat u01;// creation counts between any two given groups
  arma::mat u00;// failed creation counts between any two given groups
  arma::mat u10;// destructions counts between any two given groups
  arma::mat u11;// failed destruction counts between any two given groups
  
  // Temporary values used only during the update of an allocation
  arma::vec n_edges_to_group;// SBM part: during the update of the node, this counts the number of "successes" from a node to a group
  arma::vec n_non_edges_to_group;// SBM part: during the update of the node, this counts the number of "failures" from a node to a group
  arma::vec w01;// creation counts between a node and a group
  arma::vec w00;// failed creation counts between a node and a group
  arma::vec w10;// destructions counts between a node and a group
  arma::vec w11;// failed destruction counts between a node and a group
  
  // Hyperparameters
  double eta0, zeta0, ap, aq, bp, bq, delta;
  
  // Inference values
  double prior_value;
  arma::vec prior_value_deltas;
  double likelihood_value;
  arma::vec likelihood_value_deltas;
  double posterior_value;
  
  // GreedyICL values
  arma::vec greedy_icl_store;
  unsigned int max_n_iter;
  
  // Debug
  bool verbose;
  
  // Initialisation
  void Print();
  void Summary();
  void EvaluateActive();
  void EvaluateCountsN();
  void EvaluateNonEmptyGroups();
  void EvaluateCountsPi();
  void EvaluateCountsSBM();
  void EvaluateCountsSBTM();
  void UpdateAllValues();
  
  // Inference functions
  void EvaluatePrior();
  void EvaluateLikelihood();
  void EvaluatePosterior();
  
  // Greedy functions
  void SetUpNodeInfoForUpdate(unsigned int, unsigned int);
  void EvaluatePriorDelta(unsigned int, unsigned int, unsigned int);
  void EvaluateLikelihoodDelta(unsigned int, unsigned int, unsigned int);
  void Move(unsigned int, unsigned int, unsigned int);
  void GreedyMove(unsigned int, unsigned int);
  void GreedyOptimisation();
  void MergeUpdates();
  
  void DebugCheckAllValues();
  
protected:
private:
};

#endif // DSBTM_H
