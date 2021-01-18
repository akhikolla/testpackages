#ifndef VARSLPM_H
#define VARSLPM_H

#include <RcppArmadillo.h>
#include "utils_entropy.h"

class slpm_var
{

public:
  slpm_var(arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::cube, arma::vec, arma::vec, arma::vec, arma::vec, arma::vec, arma::vec, bool);
  void SetOptimisationPars(double, unsigned int, bool, double, double);
  void Print();
  void Summary();
  bool verbose;
  bool debug_mode;
  std::ostringstream debug_strs;

  // Global fixed quantities
  unsigned int M;
  unsigned int N;
  unsigned int K;

  // Hyperparameters
  arma::vec delta;
  double delta_sum;
  arma::vec a_gamma;
  arma::vec b_gamma;

  // Observed adjacency matrix
  arma::mat adj;

  // Variational parameters
  arma::mat var_alpha_u;
  arma::mat var_alpha_v;
  arma::mat var_beta_u;
  arma::mat var_beta_v;
  arma::cube var_lambda;
  arma::vec var_delta;
  arma::vec var_a;
  arma::vec var_b;

  // Variational statistics
  arma::vec var_lambda_sums;
  arma::vec var_s_u;
  arma::vec var_s_v;
  double var_beta_u_logs;
  double var_beta_v_logs;
  double var_delta_sum;

  // ELBO
  double term_likelihood;
  double term_prior_z;
  double term_prior_u;
  double term_prior_v;
  double term_prior_lambda;
  double term_prior_gamma;
  double term_entropy_z;
  double term_entropy_u;
  double term_entropy_v;
  double term_entropy_lambda;
  double term_entropy_gamma;
  double elbo_value;

  // Functions
  void ResetVarLambdaSums();
  void ResetVarS();
  void ResetVarBetaLogs();
  void ResetVarDeltaSum();
  void ResetTermLikelihood();
  void ResetTermPriorZ();
  void ResetTermPriorU();
  void ResetTermPriorV();
  void ResetTermPriorLambda();
  void ResetTermPriorGamma();
  void ResetTermEntropyZ();
  void ResetTermEntropyU();
  void ResetTermEntropyV();
  void ResetTermEntropyLambda();
  void ResetTermEntropyGamma();
  void ResetELBO();
  void ResetAllValues();
  void CheckValues();

  // Optimisation parameters
  arma::mat learning_rates_alpha_beta_u;
  arma::mat learning_rates_alpha_beta_v;
  double tol;
  unsigned int n_iter_max;
  bool natural_gradient;
  double learning_rate_factor_up;
  double learning_rate_factor_down;
  arma::vec elbo_values_store;

  // Optimisation functions
  arma::vec GradientU(unsigned int, unsigned int);
  arma::vec GradientV(unsigned int, unsigned int);
  void UpdateAlphaBetaU(unsigned int, unsigned int);
  void UpdateAlphaBetaV(unsigned int, unsigned int);
  void UpdateLambda(unsigned int, unsigned int);
  void UpdateA(unsigned int);
  void UpdateB(unsigned int);
  void UpdateDelta();
  void Optimisation();

protected:
private:
};

#endif // VARSLPM_H
