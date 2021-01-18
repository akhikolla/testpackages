#include "core_slpm_var.h"

arma::vec slpm_var::GradientU(unsigned int i, unsigned int k)
{
  arma::vec res;
  res.zeros(2);
  res.at(0) += - var_alpha_u.at(i,k) * var_a.at(k) / var_b.at(k);
  res.at(1) += - 0.5 * var_a.at(k) / var_b.at(k) + 0.5 / var_beta_u.at(i,k);
  for (unsigned int j=0; j<N; ++j)
  {
    double d = var_alpha_u.at(i,k) - var_alpha_v.at(j,k);
    double b = var_beta_u.at(i,k) + var_beta_v.at(j,k);
    double eta = b + d*d;
    double xi = 2*eta*eta - 2*d*d*d*d;
    res.at(0) += var_lambda.at(i,j,k) * (  R::trigamma(eta*eta/xi) * (4*d*eta/xi - 8*d*b*eta*eta/xi/xi) - 2*d/eta + 8*d*b/xi - 2*d*adj.at(i,j)  );
    res.at(1) += var_lambda.at(i,j,k) * (  R::trigamma(eta*eta/xi) * (2*eta/xi - 4*eta*eta*eta/xi/xi) - 1/eta + 4*eta/xi - adj.at(i,j)  );
  }
  return(res);
}

arma::vec slpm_var::GradientV(unsigned int j, unsigned int k)
{
  arma::vec res;
  res.zeros(2);
  res.at(0) += - var_alpha_v.at(j,k) * var_a.at(k) / var_b.at(k);
  res.at(1) += - 0.5 * var_a.at(k) / var_b.at(k) + 0.5 / var_beta_v.at(j,k);
  for (unsigned int i=0; i<M; ++i)
  {
    double d = var_alpha_v.at(j,k) - var_alpha_u.at(i,k);
    double b = var_beta_v.at(j,k) + var_beta_u.at(i,k);
    double eta = b + d*d;
    double xi = 2*eta*eta - 2*d*d*d*d;
    res.at(0) += var_lambda.at(i,j,k) * (  R::trigamma(eta*eta/xi) * (4*d*eta/xi - 8*d*b*eta*eta/xi/xi) - 2*d/eta + 8*d*b/xi - 2*d*adj.at(i,j)  );
    res.at(1) += var_lambda.at(i,j,k) * (  R::trigamma(eta*eta/xi) * (2*eta/xi - 4*eta*eta*eta/xi/xi) - 1/eta + 4*eta/xi - adj.at(i,j)  );
  }
  return(res);
}

void slpm_var::UpdateAlphaBetaU(unsigned int i, unsigned int k)
{
  double term_likelihood_delta, term_prior_u_delta, term_entropy_u_delta, elbo_value_delta;
  double alpha_old = var_alpha_u.at(i,k);
  double beta_old = var_beta_u.at(i,k);
  arma::vec gradient = GradientU(i,k);
  if (natural_gradient)
  {
    gradient.at(0) = gradient.at(0) * var_beta_u.at(i,k);
    gradient.at(1) = gradient.at(1) * 2;
  }
  
  learning_rates_alpha_beta_u.at(i,k) *= learning_rate_factor_up;
  bool stop_condition = false;
  while(!stop_condition)
  {
    double alpha_new = var_alpha_u.at(i,k) + learning_rates_alpha_beta_u.at(i,k) * gradient.at(0);
    double beta_new = var_beta_u.at(i,k) * exp(learning_rates_alpha_beta_u.at(i,k) * var_beta_u.at(i,k) * gradient.at(1));
    
    term_likelihood_delta = 0;
    for (unsigned int j=0; j<N; ++j)
    {
      double d_old = alpha_old - var_alpha_v.at(j,k);
      double d_new = alpha_new - var_alpha_v.at(j,k);
      double b_old = beta_old + var_beta_v.at(j,k);
      double b_new = beta_new + var_beta_v.at(j,k);
      double eta_old = b_old + d_old * d_old;
      double eta_new = b_new + d_new * d_new;
      double xi_old = 2 * eta_old * eta_old - 2 * d_old * d_old * d_old * d_old;
      double xi_new = 2 * eta_new * eta_new - 2 * d_new * d_new * d_new * d_new;
      term_likelihood_delta -= var_lambda.at(i,j,k) * (  R::digamma(eta_old*eta_old/xi_old) - log(eta_old) + log(xi_old) - adj.at(i,j)*eta_old  );
      term_likelihood_delta += var_lambda.at(i,j,k) * (  R::digamma(eta_new*eta_new/xi_new) - log(eta_new) + log(xi_new) - adj.at(i,j)*eta_new  );
    }
    
    term_prior_u_delta = 0;
    term_prior_u_delta -= - 0.5 * (beta_old + alpha_old*alpha_old) * var_a.at(k) / var_b.at(k);
    term_prior_u_delta += - 0.5 * (beta_new + alpha_new*alpha_new) * var_a.at(k) / var_b.at(k);
    
    term_entropy_u_delta = 0;
    term_entropy_u_delta -= 0.5 * log(2 * arma::datum::pi * beta_old);
    term_entropy_u_delta += 0.5 * log(2 * arma::datum::pi * beta_new);
    
    elbo_value_delta = term_likelihood_delta + term_prior_u_delta + term_entropy_u_delta;
    
    if (elbo_value_delta > 0)
    {
      var_alpha_u.at(i,k) = alpha_new;
      var_beta_u.at(i,k) = beta_new;
      var_s_u.at(k) -= beta_old + alpha_old * alpha_old;
      var_s_u.at(k) += beta_new + alpha_new * alpha_new;
      var_beta_u_logs -= log(2 * arma::datum::pi * beta_old);
      var_beta_u_logs += log(2 * arma::datum::pi * beta_new);
      term_likelihood += term_likelihood_delta;
      term_prior_u += term_prior_u_delta;
      term_entropy_u += term_entropy_u_delta;
      elbo_value += elbo_value_delta;
      stop_condition = true;
    }
    else
    {
      learning_rates_alpha_beta_u.at(i,k) /= learning_rate_factor_down;
    }
    if (learning_rates_alpha_beta_u.at(i,k) < 1e-6) stop_condition = true;
  }
}

void slpm_var::UpdateAlphaBetaV(unsigned int j, unsigned int k)
{
  double term_likelihood_delta, term_prior_v_delta, term_entropy_v_delta, elbo_value_delta;
  double alpha_old = var_alpha_v.at(j,k);
  double beta_old = var_beta_v.at(j,k);
  arma::vec gradient = GradientV(j,k);
  if (natural_gradient)
  {
    gradient.at(0) = gradient.at(0) * var_beta_v.at(j,k);
    gradient.at(1) = gradient.at(1) * 2;
  }
  
  learning_rates_alpha_beta_v.at(j,k) *= learning_rate_factor_up;
  bool stop_condition = false;
  while(!stop_condition)
  {
    double alpha_new = var_alpha_v.at(j,k) + learning_rates_alpha_beta_v.at(j,k) * gradient.at(0);
    double beta_new = var_beta_v.at(j,k) * exp(learning_rates_alpha_beta_v.at(j,k) * var_beta_v.at(j,k) * gradient.at(1));
    
    term_likelihood_delta = 0;
    for (unsigned int i=0; i<M; ++i)
    {
      double d_old = alpha_old - var_alpha_u.at(i,k);
      double d_new = alpha_new - var_alpha_u.at(i,k);
      double b_old = beta_old + var_beta_u.at(i,k);
      double b_new = beta_new + var_beta_u.at(i,k);
      double eta_old = b_old + d_old * d_old;
      double eta_new = b_new + d_new * d_new;
      double xi_old = 2 * eta_old * eta_old - 2 * d_old * d_old * d_old * d_old;
      double xi_new = 2 * eta_new * eta_new - 2 * d_new * d_new * d_new * d_new;
      term_likelihood_delta -= var_lambda.at(i,j,k) * (  R::digamma(eta_old*eta_old/xi_old) - log(eta_old) + log(xi_old) - adj.at(i,j)*eta_old  );
      term_likelihood_delta += var_lambda.at(i,j,k) * (  R::digamma(eta_new*eta_new/xi_new) - log(eta_new) + log(xi_new) - adj.at(i,j)*eta_new  );
    }
    
    term_prior_v_delta = 0;
    term_prior_v_delta -= - 0.5 * (beta_old + alpha_old*alpha_old) * var_a.at(k) / var_b.at(k);
    term_prior_v_delta += - 0.5 * (beta_new + alpha_new*alpha_new) * var_a.at(k) / var_b.at(k);
    
    term_entropy_v_delta = 0;
    term_entropy_v_delta -= 0.5 * log(2 * arma::datum::pi * beta_old);
    term_entropy_v_delta += 0.5 * log(2 * arma::datum::pi * beta_new);
    
    elbo_value_delta = term_likelihood_delta + term_prior_v_delta + term_entropy_v_delta;
    
    if (elbo_value_delta > 0)
    {
      learning_rates_alpha_beta_v.at(j,k) *= 2;
      var_alpha_v.at(j,k) = alpha_new;
      var_beta_v.at(j,k) = beta_new;
      var_s_v.at(k) -= beta_old + alpha_old * alpha_old;
      var_s_v.at(k) += beta_new + alpha_new * alpha_new;
      var_beta_v_logs -= log(2 * arma::datum::pi * beta_old);
      var_beta_v_logs += log(2 * arma::datum::pi * beta_new);
      term_likelihood += term_likelihood_delta;
      term_prior_v += term_prior_v_delta;
      term_entropy_v += term_entropy_v_delta;
      elbo_value += elbo_value_delta;
      stop_condition = true;
    }
    else
    {
      learning_rates_alpha_beta_v.at(j,k) /= learning_rate_factor_down;
    }
    if (learning_rates_alpha_beta_v.at(j,k) < 1e-6) stop_condition = true;
  }
}

