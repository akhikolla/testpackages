#ifndef EXPSBM_H
#define EXPSBM_H

#include <RcppArmadillo.h>

class expsbm
{

public:
  expsbm(unsigned int, arma::mat, arma::mat, arma::vec, arma::mat, arma::mat, bool, bool, double, unsigned int, bool);
  
  void Print();
  bool verbose;

  // Global fixed quantities
  unsigned int N;
  
  // Observed values stored as edgelist
  arma::mat edgelist;
  bool directed;
  bool trunc;
  
  // Observed values stored as adjacency object
  arma::field<arma::vec> A;
  arma::field<arma::vec> X;
  arma::mat W;
  void ConstructAdjacency();
  void EvaluateDataSummaries();
  
  // Clustering
  unsigned int K;
  arma::mat Z;
  arma::vec lambda;
  
  // Model parameters
  arma::mat mu, nu;
  
  // Statistics
  arma::mat A1, A0, X1, X0;// these are just summaries of the observed data -- they do not involve the model parameters -- but they are useful to make the code faster
  arma::mat L_mu, L_nu;
  arma::mat eta, zeta;
  double elbo_value;
  void EvaluateStatistics();
  void EvaluateELBO();
  
  // Optimisation
  double tol;
  unsigned int n_iter_max;
  arma::vec elbo_values_store;
  void UpdateZ(unsigned int);
  void UpdateLambda();
  void UpdateMu(unsigned int, unsigned int);
  void UpdateNu(unsigned int, unsigned int);
  void Optimisation();

protected:
private:
};

#endif // EXPSBM_H
