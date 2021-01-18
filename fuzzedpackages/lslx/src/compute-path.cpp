// Cpp class lslxOptimizer for minimizing PL criterion
// written by Po-Hsien Huang psyphh@gmail.com

#include "lslxOptimizer.h"

// [[Rcpp::depends(RcppEigen)]]
// compute solution path
// [[Rcpp::export]]
void compute_regularized_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  Rcpp::List lambda_grid = Rcpp::as<Rcpp::List>(control["lambda_grid"]);
  Rcpp::List delta_grid = Rcpp::as<Rcpp::List>(control["delta_grid"]);
  
  Rcpp::NumericVector lambda_1st_grid = Rcpp::as<Rcpp::NumericVector>(lambda_grid["lambda_1st_grid"]);
  Rcpp::NumericVector lambda_2nd_grid = Rcpp::as<Rcpp::NumericVector>(lambda_grid["lambda_2nd_grid"]);
  Rcpp::NumericVector delta_1st_grid = Rcpp::as<Rcpp::NumericVector>(delta_grid["delta_1st_grid"]);
  Rcpp::NumericVector delta_2nd_grid = Rcpp::as<Rcpp::NumericVector>(delta_grid["delta_2nd_grid"]);
  
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  int i, j, k, l, idx;
  idx = 0;
  for (i = 0; i < lambda_2nd_grid.size(); i++) {
    for (j = 0; j < delta_2nd_grid.size(); j++) {
      for (k = 0; k < lambda_1st_grid.size(); k++) {
        for (l = 0; l < delta_1st_grid.size(); l++) {
          if (!optimizer.warm_start) {
            optimizer.set_theta_value(theta_start_zero);
          }
          optimizer.set_regularizer(
            Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 
            lambda_1st_grid[k], lambda_2nd_grid[i], delta_1st_grid[l], delta_2nd_grid[j]);
          optimizer.complete_estimation();
          coefficient[idx] = optimizer.extract_coefficient();
          numerical_condition[idx] = optimizer.extract_numerical_condition();
          information_criterion[idx] = optimizer.extract_information_criterion();
          fit_index[idx] = optimizer.extract_fit_index();
          idx = idx + 1;
        }
      }
    }
  }
}



// compute stepwise solution path made by forward or backward selection
// [[Rcpp::export]]
void compute_stepwise_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  optimizer.set_regularizer(
    Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 
    0.0, 0.0, 
    INFINITY, INFINITY);
  Rcpp::IntegerVector step_grid = Rcpp::as<Rcpp::IntegerVector>(control["step_grid"]);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  int i;
  for (i = 0; i < step_grid.size(); i++) {
    if (!optimizer.warm_start) {
      optimizer.set_theta_value(theta_start_zero);
    }
    if (i == 0) {
      optimizer.complete_estimation();
      coefficient[i] = optimizer.extract_coefficient();
      numerical_condition[i] = optimizer.extract_numerical_condition();
      information_criterion[i] = optimizer.extract_information_criterion();
      fit_index[i] = optimizer.extract_fit_index();
    } else {
      optimizer.complete_searching();
      coefficient[i] = optimizer.extract_coefficient();
      numerical_condition[i] = optimizer.extract_numerical_condition();
      information_criterion[i] = optimizer.extract_information_criterion();
      fit_index[i] = optimizer.extract_fit_index();
    }
  }
}


// compute solution path made by unpenalized estimation
// [[Rcpp::export]]
void compute_none_path_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  Rcpp::NumericVector theta_start_zero = Rcpp::clone(optimizer.theta_start);
  optimizer.set_regularizer(
    Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 
    0.0, 0.0, 
    INFINITY, INFINITY);
  Rcpp::List numerical_condition = Rcpp::as<Rcpp::List>(fitted_result["numerical_condition"]);
  Rcpp::List information_criterion = Rcpp::as<Rcpp::List>(fitted_result["information_criterion"]);
  Rcpp::List fit_index = Rcpp::as<Rcpp::List>(fitted_result["fit_index"]);
  Rcpp::List coefficient = Rcpp::as<Rcpp::List>(fitted_result["coefficient"]);
  
  optimizer.complete_estimation();
  coefficient[0] = optimizer.extract_coefficient();
  numerical_condition[0] = optimizer.extract_numerical_condition();
  information_criterion[0] = optimizer.extract_information_criterion();
  fit_index[0] = optimizer.extract_fit_index();
}


// test optimization
// [[Rcpp::export]]
Rcpp::List test_optimization_cpp(
    Rcpp::List reduced_data,
    Rcpp::List reduced_model,
    Rcpp::List control,
    Rcpp::List supplied_result,
    Rcpp::List fitted_result) {
  int i;
  lslxOptimizer optimizer(reduced_data,
                          reduced_model,
                          control,
                          supplied_result);
  optimizer.set_regularizer(
    Rcpp::as< Rcpp::CharacterVector >(control["regularizer_type"]), 
    0.1, 0.0, 
    INFINITY, INFINITY);
  optimizer.update_coefficient_matrix();
  optimizer.update_implied_moment();
  optimizer.update_loss_value();
  optimizer.update_residual_weight();
  optimizer.update_model_jacobian();
  return Rcpp::wrap(optimizer.model_jacobian);
}
