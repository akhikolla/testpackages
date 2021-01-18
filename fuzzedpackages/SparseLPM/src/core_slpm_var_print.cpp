#include "core_slpm_var.h"

void slpm_var::Print()
{
  debug_strs << "\n\nclass slpm_var\n";
  debug_strs << "\nM\t=\t" << M << "\n";
  debug_strs << "\nN\t=\t" << N << "\n";
  debug_strs << "\nK\t=\t" << K << "\n";
  debug_strs << "\nAdjacency matrix:\n";
  adj.print(debug_strs);
  debug_strs << "\nvar_alpha_u:\n";
  var_alpha_u.print(debug_strs);
  debug_strs << "\nvar_alpha_v:\n";
  var_alpha_v.print(debug_strs);
  debug_strs << "\nvar_beta_u:\n";
  var_beta_u.print(debug_strs);
  debug_strs << "\nvar_beta_v:\n";
  var_beta_v.print(debug_strs);
  debug_strs << "\nvar_lambda:\n";
  var_lambda.print(debug_strs);
  debug_strs << "\nvar_delta:\n";
  var_delta.t().print(debug_strs);
  debug_strs << "\nvar_a:\n";
  var_a.t().print(debug_strs);
  debug_strs << "\nvar_b:\n";
  var_b.t().print(debug_strs);
  debug_strs << "\nStatistics:";
  debug_strs << "\nvar_lambda_sums:\n";
  var_lambda_sums.t().print(debug_strs);
  debug_strs << "\nvar_s_u:\n";
  var_s_u.t().print(debug_strs);
  debug_strs << "\nvar_s_v:\n";
  var_s_v.t().print(debug_strs);
  debug_strs << "\nvar_beta_u_logs\t=\t" << var_beta_u_logs << "\n";
  debug_strs << "\nvar_beta_v_logs\t=\t" << var_beta_v_logs << "\n";
  debug_strs << "\nvar_delta_sum\t=\t" << var_delta_sum << "\n";
  debug_strs << "\nHyperparameters:";
  debug_strs << "\ndelta:\n";
  delta.t().print(debug_strs);
  debug_strs << "\na_gamma:\n";
  a_gamma.t().print(debug_strs);
  debug_strs << "\nb_gamma:\n";
  b_gamma.t().print(debug_strs);
  debug_strs << "\nterm_likelihood\t=\t" << term_likelihood << "\n";
  debug_strs << "\nterm_prior_z\t=\t" << term_prior_z << "\n";
  debug_strs << "\nterm_prior_u\t=\t" << term_prior_u << "\n";
  debug_strs << "\nterm_prior_v\t=\t" << term_prior_v << "\n";
  debug_strs << "\nterm_prior_lambda\t=\t" << term_prior_lambda << "\n";
  debug_strs << "\nterm_prior_gamma\t=\t" << term_prior_gamma << "\n";
  debug_strs << "\nterm_entropy_z\t=\t" << term_entropy_z << "\n";
  debug_strs << "\nterm_entropy_u\t=\t" << term_entropy_u << "\n";
  debug_strs << "\nterm_entropy_v\t=\t" << term_entropy_v << "\n";
  debug_strs << "\nterm_entropy_lambda\t=\t" << term_entropy_lambda << "\n";
  debug_strs << "\nterm_entropy_gamma\t=\t" << term_entropy_gamma << "\n";
  debug_strs << "\nelbo_value\t=\t" << elbo_value << "\n";
  debug_strs << "\n\nOptimisation\n";
  debug_strs << "verbose\t=\t" << verbose << "\n";
  debug_strs << "debug_mode\t=\t" << debug_mode << "\n";
  debug_strs << "\nLearning rates U:\n";
  learning_rates_alpha_beta_u.print(debug_strs);
  debug_strs << "\nLearning rates V:\n";
  learning_rates_alpha_beta_v.print(debug_strs);
  debug_strs << "\nLearning rate modifying factor upwards\t=\t" << learning_rate_factor_up << "\n";
  debug_strs << "\nLearning rate modifying factor downwards\t=\t" << learning_rate_factor_down << "\n";
  debug_strs << "\ntol\t=\t" << tol << "\n";
  debug_strs << "Maximum number of iterations\t=\t" << n_iter_max << "\n";
  Rcpp::Rcout << debug_strs.str() << std::endl << std::endl;
}

void slpm_var::Summary()
{
  std::ostringstream strs;
  strs << "\nclass slpm_var\n";
  // strs << "\nterm_likelihood\t=\t" << term_likelihood << "\n";
  // strs << "\nterm_prior_z\t=\t" << term_prior_z << "\n";
  // strs << "\nterm_prior_u\t=\t" << term_prior_u << "\n";
  // strs << "\nterm_prior_v\t=\t" << term_prior_v << "\n";
  // strs << "\nterm_prior_lambda\t=\t" << term_prior_lambda << "\n";
  // strs << "\nterm_prior_gamma\t=\t" << term_prior_gamma << "\n";
  // strs << "\nterm_entropy_z\t=\t" << term_entropy_z << "\n";
  // strs << "\nterm_entropy_u\t=\t" << term_entropy_u << "\n";
  // strs << "\nterm_entropy_v\t=\t" << term_entropy_v << "\n";
  // strs << "\nterm_entropy_lambda\t=\t" << term_entropy_lambda << "\n";
  // strs << "\nterm_entropy_gamma\t=\t" << term_entropy_gamma << "\n";
  strs << "\nelbo_value\t=\t" << elbo_value << "\n";
  Rcpp::Rcout << strs.str() << std::endl;
}

