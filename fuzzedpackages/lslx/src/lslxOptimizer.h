// Cpp class lslxOptimizer for minimizing PL criterion
// written by Po-Hsien Huang psyphh@gmail.com

#include <iostream>
#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <string.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// define lslxOptimizer
class lslxOptimizer {
public:
  std::string loss, algorithm, regularizer_type, searcher_type;
  int iter_in_max, iter_out_max, iter_other_max, iter_armijo_max, warm_out;
  double tol_in, tol_out, tol_other;
  double step_size, armijo, momentum;
  double ridge_cov, ridge_hessian;
  bool warm_start, positive_variance, armijo_rule, enforce_cd, random_update;
  double minimum_variance;
  bool response, continuous, regularizer, searcher;

  double lambda, lambda_1st, lambda_2nd, delta, delta_1st, delta_2nd, step;
  int iter_out;
  
  int n_observation;
  Rcpp::List  sample_proportion, saturated_cov, saturated_mean, saturated_threshold, saturated_moment;
  Rcpp::List  saturated_moment_acov;
  
  int n_response, n_factor, n_eta, n_moment, n_moment_1, n_moment_2, n_group, n_theta, n_threshold;
  int idx_reference;
  Rcpp::IntegerVector idx_ordered, idx_numeric, idx_sigma, idx_gamma, idx_mu, idx_diag, idx_diag_psi, idx_nondiag; 
  Rcpp::IntegerVector idx_vech, idx_tvech, idx_vech_match, idx_nd_vech, idx_nd_tvech;
  
  Rcpp::CharacterVector theta_name;
  Rcpp::LogicalVector theta_is_free, theta_is_pen, theta_is_diag;
  Rcpp::IntegerVector theta_matrix_idx, theta_group_idx;
  Rcpp::IntegerVector theta_left_idx, theta_right_idx, theta_flat_idx, theta_tflat_idx;
  Rcpp::NumericVector theta_start, theta_value, theta_direction, theta_direction_old;
  Rcpp::LogicalVector theta_is_est, theta_is_search; 
  Rcpp::IntegerVector theta_is_est_idx, theta_is_search_idx, theta_set;
  Rcpp::CharacterVector theta_penalty;
  Rcpp::NumericVector theta_weight;
  
  double baseline_loss_value;
  int baseline_degrees_of_freedom;
  
  Eigen::MatrixXd identity_y, identity_eta, identity_theta;  
  Eigen::MatrixXd identity_y2, duplication_y;
  Eigen::MatrixXd elimination_y, duplication_eta, commutation_y;
  
  Eigen::MatrixXd some_matrix;
  
  Rcpp::List alpha, beta, beta_pinv, gamma, phi, psi;
  Rcpp::List mu, sigma, sigma_inv, implied_moment;
  Rcpp::List alpha_derivative, beta_derivative, phi_derivative;
  
  Rcpp::List model_jacobian;
  Rcpp::List model_residual;
  Rcpp::List residual_weight;
  
  double loss_value;
  Eigen::MatrixXd loss_gradient;
  Eigen::MatrixXd loss_gradient_diff;
  Eigen::MatrixXd loss_expected_hessian;
  Eigen::MatrixXd loss_observed_hessian;
  Eigen::MatrixXd loss_bfgs_hessian;
  Eigen::MatrixXd loss_bfgs_hessian_inv;
  
  double regularizer_value;
  Eigen::MatrixXd regularizer_gradient;
  
  double objective_value, objective_gradient_average;
  Eigen::MatrixXd objective_gradient;
  
  double objective_gradient_abs_max, objective_hessian_convexity;
  int n_iter_out, n_nonzero_coefficient;
  double degrees_of_freedom, robust_degrees_of_freedom, scaling_factor;
  
  double aic, aic3, caic;
  double bic, abic, hbic;
  double raic, raic3, rcaic;
  double rbic, rabic, rhbic;
  double rmsea, srmr, cfi, nnfi;
  
  lslxOptimizer(Rcpp::List reduced_data,
                Rcpp::List reduced_model,
                Rcpp::List control,
                Rcpp::List supplied_result);
  
  void set_regularizer(Rcpp::CharacterVector regularizer_type_, 
                       double lambda_1st_, double lambda_2nd_, 
                       double delta_1st_, double delta_2nd_);
  void set_searcher(Rcpp::CharacterVector searcher_type_, Rcpp::LogicalVector theta_is_search_);
  void set_theta_value(Rcpp::NumericVector theta_value_);
  
  void update_coefficient_matrix();
  void update_implied_moment();
  void update_model_jacobian();
  void update_model_jacobian_nd();
  void update_model_residual();
  void update_residual_weight();
  void update_loss_value();
  void update_loss_gradient();
  void update_loss_gradient_direct();
  void update_loss_gradient_nd();
  void update_loss_expected_hessian();
  void update_loss_observed_hessian();
  void update_loss_bfgs_hessian();
  void update_regularizer_value();
  void update_regularizer_gradient();
  void update_objective_value();
  void update_objective_gradient();
  void update_theta_direction();
  void update_nuisance();
  void update_theta_value();
  void update_theta_start();
  
  void update_coefficient();
  void update_numerical_condition();
  void update_information_criterion();
  void update_fit_index();
  
  void complete_estimation();
  void complete_searching();

  Rcpp::NumericVector extract_numerical_condition();
  Rcpp::NumericVector extract_information_criterion();
  Rcpp::NumericVector extract_fit_index();
  Rcpp::NumericVector extract_coefficient();
};
