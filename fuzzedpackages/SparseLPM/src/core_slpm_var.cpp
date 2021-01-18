#include "core_slpm_var.h"

slpm_var::slpm_var(arma::mat adj_, arma::mat var_alpha_u_, arma::mat var_alpha_v_, arma::mat var_beta_u_, arma::mat var_beta_v_, arma::cube var_lambda_, arma::vec var_delta_, arma::vec var_a_, arma::vec var_b_, arma::vec delta_, arma::vec a_gamma_, arma::vec b_gamma_, bool verbose_)
{
  adj = adj_;
  var_alpha_u = var_alpha_u_;
  var_alpha_v = var_alpha_v_;
  var_beta_u = var_beta_u_;
  var_beta_v = var_beta_v_;
  var_lambda = var_lambda_;
  var_delta = var_delta_;
  var_a = var_a_;
  var_b = var_b_;
  delta = delta_;
  delta_sum = accu(delta);
  a_gamma = a_gamma_;
  b_gamma = b_gamma_;
  M = adj.n_rows;
  N = adj.n_cols;
  K = var_alpha_u.n_cols;
  verbose = verbose_;
  debug_mode = false;
  ResetAllValues();
}

void slpm_var::SetOptimisationPars(double tol_, unsigned int n_iter_max_, bool natural_gradient_, double learning_rate_factor_up_, double learning_rate_factor_down_)
{
  tol = tol_;
  n_iter_max = n_iter_max_;
  natural_gradient = natural_gradient_;
  learning_rates_alpha_beta_u.ones(M,K);
  learning_rates_alpha_beta_v.ones(N,K);
  learning_rate_factor_up = learning_rate_factor_up_;
  learning_rate_factor_down = learning_rate_factor_down_;
}

void slpm_var::ResetVarLambdaSums()
{
  var_lambda_sums.zeros(K);
  for (unsigned int i=0; i<M; ++i) for (unsigned int j=0; j<N; ++j) for (unsigned int k=0; k<K; ++k)
  {
    var_lambda_sums.at(k) += var_lambda.at(i,j,k);
  }
}

void slpm_var::ResetVarS()
{
  var_s_u.zeros(K);
  var_s_v.zeros(K);
  for (unsigned int k=0; k<K; ++k)
  {
    for (unsigned int i=0; i<M; ++i) var_s_u.at(k) += var_beta_u.at(i,k) + var_alpha_u.at(i,k) * var_alpha_u.at(i,k);
    for (unsigned int j=0; j<N; ++j) var_s_v.at(k) += var_beta_v.at(j,k) + var_alpha_v.at(j,k) * var_alpha_v.at(j,k);
  }
}

void slpm_var::ResetVarBetaLogs()
{
  var_beta_u_logs = 0;
  var_beta_v_logs = 0;
  for (unsigned int i=0; i<M; ++i) for (unsigned int k=0; k<K; ++k)
  {
    var_beta_u_logs += log(2 * arma::datum::pi * var_beta_u.at(i,k));
  }
  for (unsigned int j=0; j<N; ++j) for (unsigned int k=0; k<K; ++k)
  {
    var_beta_v_logs += log(2 * arma::datum::pi * var_beta_v.at(j,k));
  }
}

void slpm_var::ResetVarDeltaSum()
{
  var_delta_sum = accu(var_delta);
}

void slpm_var::ResetTermLikelihood()
{
  term_likelihood = 0;
  for (unsigned int i=0; i<M; ++i) for (unsigned int j=0; j<N; ++j) for (unsigned int k=0; k<K; ++k)
  {
    double d = var_alpha_u.at(i,k) - var_alpha_v.at(j,k);
    double eta = var_beta_u.at(i,k) + var_beta_v.at(j,k) + d*d;
    double xi = 2*eta*eta - 2*d*d*d*d;
    term_likelihood += var_lambda.at(i,j,k) * ( R::digamma(eta*eta/xi) - log(eta) + log(xi) - adj.at(i,j)*eta );
  }
}

void slpm_var::ResetTermPriorZ()
{
  term_prior_z = 0;
  for (unsigned int k=0; k<K; ++k)
  {
    term_prior_z += (R::digamma(var_delta.at(k)) - R::digamma(var_delta_sum)) * var_lambda_sums.at(k);
  }
}

void slpm_var::ResetTermPriorU()
{
  term_prior_u = -0.5*M*K*log(2*arma::datum::pi);
  for (unsigned int k=0; k<K; ++k)
  {
    term_prior_u += 0.5 * M * (R::digamma(var_a.at(k)) - log(var_b.at(k)));
    term_prior_u -= 0.5 * var_s_u.at(k) * var_a.at(k) / var_b.at(k);
  }
}

void slpm_var::ResetTermPriorV()
{
  term_prior_v = -0.5*N*K*log(2*arma::datum::pi);
  for (unsigned int k=0; k<K; ++k)
  {
    term_prior_v += 0.5 * N * (R::digamma(var_a.at(k)) - log(var_b.at(k)));
    term_prior_v -= 0.5 * var_s_v.at(k) * var_a.at(k) / var_b.at(k);
  }
}

void slpm_var::ResetTermPriorLambda()
{
  term_prior_lambda = lgamma(delta_sum);
  for (unsigned int k=0; k<K; ++k)
  {
    term_prior_lambda -= lgamma(delta.at(k));
    term_prior_lambda += (delta.at(k)-1) * (R::digamma(var_delta.at(k)) - R::digamma(var_delta_sum));
  }
}

void slpm_var::ResetTermPriorGamma()
{
  term_prior_gamma = 0;
  for (unsigned int k=0; k<K; ++k)
  {
    term_prior_gamma += a_gamma.at(k) * log(b_gamma.at(k)) - lgamma(a_gamma.at(k));
    term_prior_gamma += (a_gamma.at(k)-1) * (R::digamma(var_a.at(k)) - log(var_b.at(k)));
    term_prior_gamma += - b_gamma.at(k) * var_a.at(k) / var_b.at(k);
  }
}

void slpm_var::ResetTermEntropyZ()
{
  term_entropy_z = 0;
  for (unsigned int i=0; i<M; ++i) for (unsigned int j=0; j<N; ++j) for (unsigned int k=0; k<K; ++k)
  {
    term_entropy_z += Entropy(var_lambda.at(i,j,k));
  }
}

void slpm_var::ResetTermEntropyU()
{
  term_entropy_u = 0.5*M*K + 0.5 * var_beta_u_logs;
}

void slpm_var::ResetTermEntropyV()
{
  term_entropy_v = 0.5*N*K + 0.5 * var_beta_v_logs;
}

void slpm_var::ResetTermEntropyLambda()
{
  term_entropy_lambda = - lgamma(var_delta_sum);
  for (unsigned int k=0; k<K; ++k)
  {
    term_entropy_lambda += lgamma(var_delta.at(k));
    term_entropy_lambda -= (var_delta.at(k) - 1) * (R::digamma(var_delta.at(k)) - R::digamma(var_delta_sum));
  }
}

void slpm_var::ResetTermEntropyGamma()
{
  term_entropy_gamma = 0;
  for (unsigned int k=0; k<K; ++k)
  {
    term_entropy_gamma += - var_a.at(k) * R::digamma(var_a.at(k)) + R::digamma(var_a.at(k));
    term_entropy_gamma += var_a.at(k) + lgamma(var_a.at(k)) - log(var_b.at(k));
  }
}

void slpm_var::ResetELBO()
{
  elbo_value = 0;
  elbo_value += term_likelihood;
  elbo_value += term_prior_z + term_prior_u + term_prior_v + term_prior_lambda + term_prior_gamma;
  elbo_value += term_entropy_z + term_entropy_u + term_entropy_v + term_entropy_lambda + term_entropy_gamma;
}

void slpm_var::ResetAllValues()
{
  ResetVarLambdaSums();
  ResetVarS();
  ResetVarBetaLogs();
  ResetVarDeltaSum();
  ResetTermLikelihood();
  ResetTermPriorZ();
  ResetTermPriorLambda();
  ResetTermPriorU();
  ResetTermPriorV();
  ResetTermPriorGamma();
  ResetTermEntropyZ();
  ResetTermEntropyU();
  ResetTermEntropyV();
  ResetTermEntropyLambda();
  ResetTermEntropyGamma();
  ResetELBO();
}
