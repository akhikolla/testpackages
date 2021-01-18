#include "core_slpm_var.h"

void slpm_var::UpdateLambda(unsigned int i, unsigned int j)
{
  arma::vec log_terms_likelihood, log_terms_prior;
  log_terms_likelihood.zeros(K);
  log_terms_prior.zeros(K);
  for (unsigned int k=0; k<K; ++k)
  {
    double d = var_alpha_u.at(i,k) - var_alpha_v.at(j,k);
    double eta = var_beta_u.at(i,k) + var_beta_v.at(j,k) + d*d;
    double xi = 2*eta*eta - 2*d*d*d*d;
    log_terms_prior.at(k) = R::digamma(var_delta.at(k)) - R::digamma(var_delta_sum);
    log_terms_likelihood.at(k) = R::digamma(eta*eta/xi) - log(eta) + log(xi) - adj.at(i,j)*eta;
  }
  
  for (unsigned int k=0; k<K; ++k) 
  {
    term_likelihood -= var_lambda.at(i,j,k) * log_terms_likelihood.at(k);
    term_prior_z -= var_lambda.at(i,j,k) * log_terms_prior.at(k);
    term_entropy_z -= Entropy(var_lambda.at(i,j,k));
    var_lambda_sums.at(k) -= var_lambda.at(i,j,k);
  }
  
  double prop_const = max(log_terms_prior + log_terms_likelihood);
  for (unsigned int k=0; k<K; ++k) var_lambda.at(i,j,k) = exp(log_terms_likelihood.at(k) + log_terms_prior.at(k) - prop_const);
  double values_sum = 0;
  for (unsigned int k=0; k<K; ++k) values_sum += var_lambda.at(i,j,k);
  for (unsigned int k=0; k<K; ++k) var_lambda.at(i,j,k) = var_lambda.at(i,j,k) / values_sum;
  
  for (unsigned int k=0; k<K; ++k) 
  {
    term_likelihood += var_lambda.at(i,j,k) * log_terms_likelihood.at(k);
    term_prior_z += var_lambda.at(i,j,k) * log_terms_prior.at(k);
    term_entropy_z += Entropy(var_lambda.at(i,j,k));
    var_lambda_sums.at(k) += var_lambda.at(i,j,k);
  }
  ResetELBO();
}

void slpm_var::UpdateA(unsigned int k)
{
  term_prior_u -= 0.5 * M * R::digamma(var_a.at(k)) - 0.5 * var_a.at(k) * var_s_u.at(k) / var_b.at(k);
  term_prior_v -= 0.5 * N * R::digamma(var_a.at(k)) - 0.5 * var_a.at(k) * var_s_v.at(k) / var_b.at(k);
  term_prior_gamma -= (a_gamma.at(k)-1) * R::digamma(var_a.at(k)) - b_gamma.at(k) * var_a.at(k) / var_b.at(k);
  term_entropy_gamma += (var_a.at(k)-1) * R::digamma(var_a.at(k)) - var_a.at(k) - lgamma(var_a.at(k));
  
  var_a.at(k) = a_gamma.at(k) + 0.5 * (M + N);
  
  term_prior_u += 0.5 * M * R::digamma(var_a.at(k)) - 0.5 * var_a.at(k) * var_s_u.at(k) / var_b.at(k);
  term_prior_v += 0.5 * N * R::digamma(var_a.at(k)) - 0.5 * var_a.at(k) * var_s_v.at(k) / var_b.at(k);
  term_prior_gamma += (a_gamma.at(k)-1) * R::digamma(var_a.at(k)) - b_gamma.at(k) * var_a.at(k) / var_b.at(k);
  term_entropy_gamma -= (var_a.at(k)-1) * R::digamma(var_a.at(k)) - var_a.at(k) - lgamma(var_a.at(k));
  ResetELBO();
}

void slpm_var::UpdateB(unsigned int k)
{
  term_prior_u -= - 0.5 * M * log(var_b.at(k)) - 0.5 * var_a.at(k) * var_s_u.at(k) / var_b.at(k);
  term_prior_v -= - 0.5 * N * log(var_b.at(k)) - 0.5 * var_a.at(k) * var_s_v.at(k) / var_b.at(k);
  term_prior_gamma -= - (a_gamma.at(k)-1) * log(var_b.at(k)) - b_gamma.at(k) * var_a.at(k) / var_b.at(k);
  term_entropy_gamma -= - log(var_b.at(k));
  
  var_b.at(k) = b_gamma.at(k) + 0.5 * var_s_u.at(k) + 0.5 * var_s_v.at(k);
  
  term_prior_u += - 0.5 * M * log(var_b.at(k)) - 0.5 * var_a.at(k) * var_s_u.at(k) / var_b.at(k);
  term_prior_v += - 0.5 * N * log(var_b.at(k)) - 0.5 * var_a.at(k) * var_s_v.at(k) / var_b.at(k);
  term_prior_gamma += - (a_gamma.at(k)-1) * log(var_b.at(k)) - b_gamma.at(k) * var_a.at(k) / var_b.at(k);
  term_entropy_gamma += - log(var_b.at(k));
  ResetELBO();
}

void slpm_var::UpdateDelta()
{
  for (unsigned int k=0; k<K; ++k) var_delta.at(k) = delta.at(k) + var_lambda_sums.at(k);
  ResetVarDeltaSum();
  ResetTermPriorZ();
  ResetTermPriorLambda();
  ResetTermEntropyLambda();
  ResetELBO();
}



