#include "core_slpm_var.h"

void slpm_var::CheckValues()
{
  double term_likelihood_check = term_likelihood;
  double term_prior_z_check = term_prior_z;
  double term_prior_u_check = term_prior_u;
  double term_prior_v_check = term_prior_v;
  double term_prior_lambda_check = term_prior_lambda;
  double term_prior_gamma_check = term_prior_gamma;
  double term_entropy_z_check = term_entropy_z;
  double term_entropy_u_check = term_entropy_u;
  double term_entropy_v_check = term_entropy_v;
  double term_entropy_lambda_check = term_entropy_lambda;
  double term_entropy_gamma_check = term_entropy_gamma;
  double elbo_value_check = elbo_value;
  ResetAllValues();
  Rcpp::Rcout << "\n\n----------------- DEBUG START -----------------" << std::endl;
  Rcpp::Rcout << "error on term_likelihood\t = \t" << std::abs(term_likelihood_check - term_likelihood) << std::endl;
  Rcpp::Rcout << "error on term_prior_z\t\t = \t" << std::abs(term_prior_z_check - term_prior_z) << std::endl;
  Rcpp::Rcout << "error on term_prior_u\t\t = \t" << std::abs(term_prior_u_check - term_prior_u) << std::endl;
  Rcpp::Rcout << "error on term_prior_v\t\t = \t" << std::abs(term_prior_v_check - term_prior_v) << std::endl;
  Rcpp::Rcout << "error on term_prior_lambda\t = \t" << std::abs(term_prior_lambda_check - term_prior_lambda) << std::endl;
  Rcpp::Rcout << "error on term_prior_gamma\t = \t" << std::abs(term_prior_gamma_check - term_prior_gamma) << std::endl;
  Rcpp::Rcout << "error on term_entropy_z\t\t = \t" << std::abs(term_entropy_z_check - term_entropy_z) << std::endl;
  Rcpp::Rcout << "error on term_entropy_u\t\t = \t" << std::abs(term_entropy_u_check - term_entropy_u) << std::endl;
  Rcpp::Rcout << "error on term_entropy_v\t\t = \t" << std::abs(term_entropy_v_check - term_entropy_v) << std::endl;
  Rcpp::Rcout << "error on term_entropy_lambda\t = \t" << std::abs(term_entropy_lambda_check - term_entropy_lambda) << std::endl;
  Rcpp::Rcout << "error on term_entropy_gamma\t = \t" << std::abs(term_entropy_gamma_check - term_entropy_gamma) << std::endl;
  Rcpp::Rcout << "error on elbo_value\t\t = \t" << std::abs(elbo_value_check - elbo_value) << std::endl;
  Rcpp::Rcout << "----------------- DEBUG END -------------------\n\n" << std::endl;
}

