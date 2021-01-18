// Cpp class lslxOptimizer for minimizing PL criterion
// written by Po-Hsien Huang psyphh@gmail.com

#include "lslxOptimizer.h"
#include "utility-function.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// compute coefficient matrix
// [[Rcpp::export]]
Rcpp::List compute_coefficient_matrix_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List coefficient_matrix;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  
  coefficient_matrix = 
    Rcpp::List::create(Rcpp::Named("alpha") = optimizer.alpha,
                       Rcpp::Named("beta") = optimizer.beta,
                       Rcpp::Named("phi") = optimizer.phi);
  
  return Rcpp::wrap(coefficient_matrix);
}

// compute implied covariance
// [[Rcpp::export]]
Rcpp::List compute_implied_cov_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List implied_cov;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  
  implied_cov = optimizer.sigma;
  return Rcpp::wrap(implied_cov);
}

// compute implied mean
// [[Rcpp::export]]
Rcpp::List compute_implied_mean_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Rcpp::List implied_mean;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  
  implied_mean = optimizer.mu;
  return Rcpp::wrap(implied_mean);
}

// compute moment jacobian
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_model_jacobian_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Eigen::MatrixXd model_jacobian = 
    Eigen::MatrixXd::Zero(optimizer.n_group * optimizer.n_moment, 
                          optimizer.n_theta);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  optimizer.update_model_jacobian();
  int i;
  for (i = 0; i < optimizer.n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> model_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.model_jacobian[i]));
    model_jacobian.block(
      i * optimizer.n_moment, 0,
      optimizer.n_moment, optimizer.n_theta) = model_jacobian_i;
  }
  return Rcpp::wrap(model_jacobian);
}



// compute expected fisher
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_expected_information_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd expected_information;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  
  optimizer.update_residual_weight();
  optimizer.update_model_jacobian();
  optimizer.update_loss_expected_hessian();
  expected_information = 0.5 * optimizer.loss_expected_hessian;
  return Rcpp::wrap(expected_information);
}

// compute observed fisher
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_observed_information_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd observed_information;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_theta_start();
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  
  optimizer.update_residual_weight();
  optimizer.update_model_jacobian();
  optimizer.update_loss_gradient();
  optimizer.update_loss_observed_hessian();
  observed_information = 0.5 * optimizer.loss_observed_hessian;
  return Rcpp::wrap(observed_information);
}

// compute asymptotic covariance of score
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_score_acov_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Eigen::MatrixXd score_acov = 
    Eigen::MatrixXd::Zero(optimizer.n_theta, 
                          optimizer.n_theta);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  optimizer.update_residual_weight();
  optimizer.update_model_jacobian();
  
  int i;
  for (i = 0; i < optimizer.n_group; i++) {
    Eigen::Map<Eigen::MatrixXd> model_jacobian_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.model_jacobian[i]));
    Eigen::Map<Eigen::MatrixXd> residual_weight_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.residual_weight[i]));
    Eigen::Map<Eigen::MatrixXd> saturated_moment_acov_i(Rcpp::as< Eigen::Map <Eigen::MatrixXd> >(optimizer.saturated_moment_acov[i]));
    score_acov += 
      model_jacobian_i.transpose() * residual_weight_i * 
      saturated_moment_acov_i * residual_weight_i * model_jacobian_i;
  }
  return Rcpp::wrap(score_acov);
}

// compute loss value
// [[Rcpp::export]]
double compute_loss_value_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  double loss_value;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  optimizer.update_loss_value();
  loss_value = optimizer.loss_value;
  return loss_value;
}

// compute loss gradient 
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_loss_gradient_cpp(
    Rcpp::NumericVector theta_value,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd loss_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  if (optimizer.loss == "ml") {
    optimizer.update_loss_gradient_direct();
  } else {    
    optimizer.update_residual_weight();
    optimizer.update_model_jacobian();
    optimizer.update_loss_gradient();
    }
  loss_gradient = optimizer.loss_gradient;
  return Rcpp::wrap(loss_gradient);
}

// compute regularizer gradient
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_regularizer_gradient_cpp(
    Rcpp::NumericVector theta_value,
    double lambda_1st,
    double lambda_2nd,
    double delta_1st,
    double delta_2nd,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd regularizer_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<Rcpp::CharacterVector>(control["regularizer_type"]), 
                            lambda_1st, lambda_2nd, delta_1st, delta_2nd);
  optimizer.update_regularizer_gradient();
  regularizer_gradient = optimizer.regularizer_gradient;
  return Rcpp::wrap(regularizer_gradient);
}

// compute objective gradient
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_objective_gradient_cpp(
    Rcpp::NumericVector theta_value,
    double lambda_1st,
    double lambda_2nd,
    double delta_1st,
    double delta_2nd,
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result) {
  Eigen::MatrixXd objective_gradient;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_theta_value(theta_value);
  optimizer.set_regularizer(Rcpp::as<Rcpp::CharacterVector>(control["regularizer_type"]), 
                            lambda_1st, lambda_2nd, delta_1st, delta_2nd);
  
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  
  if (optimizer.loss == "ml") {
    optimizer.update_loss_gradient_direct();
  } else {    
    optimizer.update_residual_weight();
    optimizer.update_model_jacobian();
    optimizer.update_loss_gradient();
  }
  optimizer.update_regularizer_gradient();
  optimizer.update_objective_gradient();
  objective_gradient = optimizer.objective_gradient;
  return Rcpp::wrap(objective_gradient);
}
