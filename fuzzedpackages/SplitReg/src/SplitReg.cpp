#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
#include "config.h"


// Function to return the positive part of any number
double Positive_Part(const double & x){
  
  double positive = x;
  if (x <= 0){
    positive = 0; 
  }
  return positive;
}

// Function to return the sign of any number - returns a sign of 0 if the numerical argument is 0
double Sign(const double & x){
  
  double sign = 0;
  if (x < 0){
    sign = -1; 
  } 
  if (x > 0){
    sign = 1;
  }
  return sign;
}

// Function that returns the absolute value of any number
double Absolute_Value(const double & x){
  
  double abs_value = x * Sign(x);
  return abs_value;
}

// Function that returns the numerical value from the soft-thresholding operator (takes 2 numerical values as input)
double Soft_Thresholding(const double & x, 
                         const double & gamma){
  
  double soft = 0;
  soft = Sign(x) * Positive_Part(Absolute_Value(x) - gamma);
  return soft;
}


arma::mat beta_weights(const arma::mat & beta,
                       const arma::uword & group){
  // Computes weights for the l1 interaction penalty term
  arma::uword num_groups = beta.n_cols;
  arma::mat sum_abs = zeros(beta.n_rows, 1);
  arma::vec indices = ones(num_groups, 1);
  indices[group] = 0;
  sum_abs = abs(beta) * indices;
  return(sum_abs);
}

double EN_penalty(const arma::mat & beta,
                  const double & lambda_sparsity,
                  const double & alpha){
  // Compute the EN penalty
  double penalty = lambda_sparsity * ((1 - alpha) * 0.5 * pow(norm(beta, "fro"), 2) + alpha * accu(abs(beta)));
  return(penalty);
}

double Diversity_Penalty(const arma::mat & beta,
                         const double & lambda_diversity){
  // Function to compute the diversity penalty
  double penalty = 0;
  arma::mat gram_beta = zeros(beta.n_rows, beta.n_rows);
  gram_beta = abs(beta.t()) * abs(beta);
  gram_beta.diag().zeros();
  penalty = 0.5 * accu(gram_beta);
  penalty *= lambda_diversity;
  return(penalty);
}

// [[Rcpp::export]]
double Ensemble_EN_Objective(const arma::mat & current_res,
                             const arma::mat & beta,
                             const double & lambda_sparsity,
                             const double & lambda_diversity,
                             const double & alpha){
  // Compute the Ensemble EN objective function
  double n = current_res.n_rows;
  arma::mat squared_res = square(current_res);
  double loss = accu(squared_res / (2 * n));
  double EN_pen = EN_penalty(beta, lambda_sparsity, alpha);
  double ensemble_EN_pen = Diversity_Penalty(beta, lambda_diversity);
  double objective = loss + EN_pen + ensemble_EN_pen;
  return objective;
}

// [[Rcpp::export]]
arma::cube Prediction_Grid(const arma::mat & x_test,
                           const arma::mat & x_train,
                           const arma::vec & y_train,
                           const arma::cube & grid_betas){
  // Function that returns predictions from a sequence of betas (coefficients)
  arma::uword n = x_test.n_rows;
  arma::uword len_grid = grid_betas.n_slices;
  arma::uword num_models = grid_betas.n_cols;
  arma::cube predictions = zeros(n, num_models, len_grid);
  arma::rowvec mu_x = mean(x_train);
  double mu_y = mean(y_train);
  for(arma::uword i = 0; i < len_grid; i++){
    predictions.slice(i) =  mu_y + x_test * grid_betas.slice(i);
    predictions.slice(i).each_row() -= mu_x * grid_betas.slice(i);
  }
  return(predictions);
}

void Cycling(const arma::mat & x,
             const arma::vec & y,
             const arma::mat & thresh,
             const double & stdz,
             const arma::uword & group,
             arma::mat & current_res,
             arma::mat & out_beta){
  
  // Does one cycle of CD for one group, updates the coefficient
  // in out_beta and the residuals in current_res
  double n = x.n_rows;
  arma::uword p = x.n_cols;
  double resid_corr = 0;
  double old_coef = 0;
  
  for(arma::uword j = 0; j < p; j++){
    old_coef = out_beta(j, group);
    out_beta(j, group) = 0;
    // Current residuals
    resid_corr = dot(x.col(j), current_res.col(group) / n) + old_coef;
    // Update
    out_beta(j, group) = Soft_Thresholding(resid_corr, thresh[j]) / stdz;
    if (out_beta(j, group) != old_coef){
      // Only update residuals if the coefficient changed
      current_res.col(group) += x.col(j) * (old_coef - out_beta(j, group));
    } 
  }
}


void Ensemble_EN_Solver(const arma::mat & x,
                        const arma::vec & y,
                        const double & lambda_sparsity,
                        const double & lambda_diversity,
                        const double & alpha,
                        const arma::uword & num_groups,
                        const double & tolerance,
                        unsigned long & max_iter,
                        arma::mat & current_res,
                        arma::mat & beta){
  // Solves ensembles EN function for fixed penalty terms. Assumes x and y
  // have been centered / scaled. output is not de-standardized
  // Input
  // x: design matrix, centered to zero mean and scaled to unit variance
  // y: responses, centered to zero mean and scale to unit variance
  // lambdas_sparsity: penalty parameter for individual coefficients
  // lambdas_diversity: penalty parameter for interactions between groups
  // alpha: Elastic Net tuning constant
  // num_groups: number of groups
  // tolerance: tolerance parameter to stop the iterations
  // max_iter: maximum number of iterations before stopping the iterations over the groups
  //   
  // # Output
  // beta: slopes
  // current_res: residuals
  arma::uword p = x.n_cols;
  arma::mat thresh = zeros(p, 1);
  double stdz = 0;
  double conv_crit = 1;
  arma::uword iteration = 0;
  arma::mat beta_old = zeros(p, num_groups);
  
  beta_old = beta;
  stdz = 1 + lambda_sparsity * (1 - alpha);
  // Do one cycle to start with
  iteration += 1;
  for (arma::uword group = 0; group < num_groups; group++){
    // Update penalty
    thresh = lambda_sparsity * alpha + lambda_diversity * beta_weights(beta, group);
    // Do one CD cycle
    Cycling(x, y, thresh, stdz, group, current_res, beta);
  }
  beta_old = beta;
  // cout << '\n' << objective_new << '\n';
  while((conv_crit > tolerance) & (iteration <= max_iter)){
    iteration += 1;
     for (arma::uword group = 0; group < num_groups; group++){
       // Update penalty
       thresh = lambda_sparsity * alpha + lambda_diversity * beta_weights(beta, group);
       // Do one CD cycle
       Cycling(x, y, thresh, stdz, group, current_res, beta);
    }
    conv_crit = square(mean(beta_old, 1) - mean(beta, 1)).max();
    // objective_old = objective_new;
    beta_old = beta;
  }
}

// [[Rcpp::export]]
arma::cube Ensemble_EN_Grid(const arma::mat & x,
                            const arma::vec & y,
                            const int & which_lambda,
                            const arma::vec & lambdas_grid,
                            const double & lambda_fixed,
                            const double & alpha,
                            const arma::uword & num_groups,
                            const double & tolerance,
                            unsigned long & max_iter){
  // Computes Ensemble EN over a path of penalty values
  //   
  // Input
  // x: design matrix
  // y: responses
  // which_lambda: which penalty is the grid for? 1: lambda_sparsity, 2: lambda_diversity
  // lambdas_grid: grid of penalty values to compute the solution over
  // lambda_fixed: the other penalty
  // alpha: EN constant
  // num_groups: number of groups
  // tolerance: tolerance parameter to stop the iterations
  // max_iter: maximum number of iterations before stopping
  //   
  // Output
  // a cube whose slices are the slopes computer over lambda_grid
  
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  arma::mat x_std = x;
  arma::vec y_std = y;
  arma::mat mu_x = zeros(1, p);
  arma::mat sd_x = zeros(1, p);
  arma::uword num_lambda = lambdas_grid.n_elem;
  double sd_y = 0;
  double mu_y = 0;
  // Slopes
  arma::mat beta_old_grid = zeros(p, num_groups);
  // Current model residuals for each group (by column), at each iteration and grid point
  arma::mat current_res = zeros(n, num_groups);
  // Output
  arma::cube out_betas = zeros(p, num_groups, num_lambda);
  mu_x = mean(x_std);
  sd_x = stddev(x_std, 1);
  x_std.each_row() -= mu_x;
  x_std.each_row() /= sd_x;
  mu_y = mean(y_std);
  y_std = y - mu_y;
  sd_y = stddev(y_std, 1);
  y_std = y / sd_y;
  // Residuals are centered and scaled reponses for the empty model
  current_res.each_col() = y_std;
  // sd_x is not the sd of x
  sd_x = sd_x / sd_y;
  
  if (which_lambda==1){
    for (int i = (num_lambda - 1); i >= 0; i--){
      // Use the solver. Iterations start at beta_old_grid. Output
      // is written to beta_old_grid, residuals are updated in current_res
      Ensemble_EN_Solver(x_std, y_std, lambdas_grid[i], lambda_fixed,
                        alpha, num_groups, tolerance, max_iter, current_res, beta_old_grid);
      out_betas.slice(i) = beta_old_grid;
      // De-standardization of the beta coefficients
      out_betas.slice(i).each_col() /= (sd_x.t());
    }
  } else{
      for (int i = (num_lambda - 1); i >= 0; i--){
        // Use the solver. Iterations start at beta_old_grid. 
        // Output is written to beta_old_grid
        Ensemble_EN_Solver(x_std, y_std, lambda_fixed, lambdas_grid[i],
                          alpha, num_groups, tolerance, max_iter, current_res, beta_old_grid);
        out_betas.slice(i) = beta_old_grid;
        // De-standardization f the beta coefficients
        out_betas.slice(i).each_col() /= (sd_x.t());
    }
  }
  
  return(out_betas);
}

bool Check_Interactions_Beta(const arma::mat & beta){
  // This function checks if there are interactions between groups in the matrix of betas
  arma::uword p = beta.n_rows;
  bool interactions = false;
  for (arma::uword i = 0; i < p; i ++){
    arma::mat temp = nonzeros(beta.row(i));
    if (temp.n_rows > 1){
      return(true);
    }
  }
  return(interactions);
}

arma::uvec Check_Interactions(const arma::cube & betas){
  // Returns a vector with ones corresponding to the betas that have interactions.
  arma::vec checks = zeros(betas.n_slices, 1);
  arma::vec all_ones = ones(betas.n_slices, 1);
  for(arma::uword i = 0; i < betas.n_slices; i++){
    checks(i) = Check_Interactions_Beta(betas.slice(i));
  }
  return(checks==all_ones);
}

arma::vec Lambdas_Diversity_Grid(const arma::mat & x,
                                 const arma::vec & y,
                                 const double & lambda_sparsity_min,
                                 const arma::uword & len_grid,
                                 const double & alpha,
                                 const double & eps,
                                 const arma::uword & num_groups,
                                 const double & tolerance,
                                 unsigned long & max_iter){
  // Finds a reasonable grid of lambda_diversity by attempting to find the smallest lambda_diversity that kills 
  // all interactions at lambda_sparsity_min.
  // May be improved in the future.   
  //
  // Input
  // x: design matrix, centered and scaled
  // y: responses, centered and scaled
  // lambdas_sparsity_min: smallest lambda_in to be used
  // len_grid: length of the desired grid
  // alpha: EN tuning constant
  // eps: ratio of lambda_diversity_min to lambda_diversity max
  // num_groups: number of groups
  // tolerance: tolerance parameter to stop the iterations
  // max_iter: maximum number of iterations before stopping
  //   
  // Output
  // lambdas_grid: a grid of lambda ints with length len_grid
  
  double p = x.n_cols;
  arma::cube betas = zeros(p, num_groups, len_grid);
  // Initial guess for the diversity penalty
  double lambdas_diversity_max = num_groups;
  arma::vec lambdas_grid = exp(linspace(log(eps * lambdas_diversity_max), log(lambdas_diversity_max), len_grid));
  betas = Ensemble_EN_Grid(x, y, 2, lambdas_grid, lambda_sparsity_min, alpha, num_groups, tolerance, max_iter);
  arma::uword counter = 0;
  // While interactions remain, increase lambdas_diversity_max by scaling it by a constant factor of two
  while (Check_Interactions_Beta(betas.slice(len_grid - 1)) & (counter <= 5)){
    counter += 1;
    lambdas_diversity_max = lambdas_diversity_max * 2;
    lambdas_grid = exp(linspace(log(eps * lambdas_diversity_max), log(lambdas_diversity_max), len_grid));
    betas = Ensemble_EN_Grid(x, y, 2, lambdas_grid, lambda_sparsity_min, alpha, num_groups, tolerance, max_iter);
  }
  // If we could not kill all the interactions
  if (Check_Interactions_Beta(betas.slice(len_grid - 1))){
    Rcpp::warning("Failure to find lambda_diversity that kills all interactions");
    return(lambdas_grid);
  } else {
    // Find smallest lambda_diversity in the grid such that there are no interactions
    arma::uvec interactions = Check_Interactions(betas);
    // Find smallest index where there are no interactions
    arma::uvec indexes = find(interactions==0, 1);
    arma::uword index = indexes[0];
    lambdas_diversity_max = lambdas_grid[index];
    lambdas_grid = exp(linspace(log(eps * lambdas_diversity_max), log(lambdas_diversity_max), len_grid));
  }
  lambdas_grid.insert_rows(0, 1);
  return(lambdas_grid);
}

arma::uvec Set_Diff(const arma::uvec & big,
                    const arma::uvec & small){
  // Find set difference between a big and a small set of variables.
  // Note: small is a subset of big (both are sorted).
  int m = small.n_elem;
  int n = big.n_elem;
  arma::uvec test = uvec(n, fill::zeros);
  arma::uvec zeros = uvec(n - m, fill::zeros);
  
  for (int j = 0 ; j < m ; j++){
    test[small[j]] = small[j];
  }
  
  test = big - test;
  if(small[0] != 0){
    test[0] = 1;
  }
  zeros = find(test != 0);
  return(zeros);
}
// [[Rcpp::export]]
arma::vec CV_Ensemble_EN(const arma::mat & x,
                         const arma::vec & y,
                         const arma::uword & which_lambda,
                         const arma::vec & lambdas_grid,
                         const double & lambda_fixed,
                         const double & alpha,
                         const arma::uword & num_groups,
                         const arma::uword & num_folds,
                         const double & tolerance,
                         unsigned long & max_iter,
                         const arma::uword & num_threads){
  // Finds CV MSE for Ensemble EN with given penalty parameters
  // Input
  // x: design matrix, shuffled
  // y: responses, shuffled
  // which_lambda: which penalty is the grid for? 1: lambda_sparsity, 2: lambda_diversity
  // lambdas_grid: grid of penalty values to compute the solution over
  // lambda_fixed: the other penalty
  // alpha: EN tuning
  // num_groups: number of groups
  // num_folds: number of folds for CV
  // tolerance: tolerance parameter to stop the iterations
  // max_iter: maximum number of iterations before stopping
  // num_threads: number of threads for parallel computations
  //   
  // Output
  // mse : the CV MSE for each lambda in lambda_grid
  const arma::uword p = x.n_cols;
  const double n = x.n_rows;
  const arma::uword num_lambdas = lambdas_grid.n_rows;
  const arma::uvec indin = linspace<uvec>(0, n - 1, n);
  const arma::uvec inint = linspace<uvec>(0, n , num_folds + 1);
  arma::mat mses = zeros(num_lambdas, num_folds);
# pragma omp parallel for num_threads(num_threads)
  for(arma::uword fold = 0; fold < num_folds; fold++){
    // Get test and training samples
    arma::uvec test = linspace<uvec>(inint[fold], inint[fold + 1] - 1, inint[fold + 1] - inint[fold]);
    arma::uvec train = Set_Diff(indin, test);
    // Fit using train, predict using test
    arma::cube betas = zeros(p, num_groups, num_lambdas);
    betas = Ensemble_EN_Grid(x.rows(train), y.rows(train), which_lambda,
                             lambdas_grid, lambda_fixed, alpha, num_groups,
                             tolerance, max_iter);
    arma::cube preds = Prediction_Grid(x.rows(test), x.rows(train), y.rows(train), betas);
    arma::mat preds_ave = mean(preds, 1);
    for(arma::uword i = 0; i < num_lambdas; i++){
      mses.at(i, fold) = accu(square(y.rows(test)/sqrt(n) - preds_ave.col(i)/sqrt(n)));
    }
  }
  arma::vec out = sum(mses, 1);
  return(out);
}

// [[Rcpp::export]]
List Main_Ensemble_EN(const arma::mat & x_perm,
                      const arma::vec & y_perm,
                      const arma::uword num_lambdas_sparsity,
                      const arma::uword num_lambdas_diversity,
                      const double & alpha,
                      const arma::uword & num_groups,
                      const double & tolerance,
                      unsigned long & max_iter,
                      const arma::uword & num_folds,
                      const arma::uword & num_threads){
  // Finds optimal penalties for a Ensemble EN in a sequential fashion. 
  // Optimality is measure by CV MSE.
  //   
  // Input
  // x_perm: design matrix, rows shuffled
  // y_perm: responses, rows shuffled
  // num_lambdas_sparsity: number of penalty parameters for individual coefficients
  // num_lambdas_diversity: number of penalty parameters for interactions between groups
  // alpha: EN tuning constant
  // num_groups: number of groups
  // tolerance: tolerance parameter to stop the iterations
  // max_iter: maximum number of iterations before stopping
  // num_folds: number of folds for CV
  // num_threads: number of threads for parallel computations
  //  
  // Output
  // a list with entries
  // betas: betas computed over lambdas_sparsity whole path, with lambda_diversity fixed at the optimum
  // index_opt: index of the optimal lambda_sparsity
  // lambda_sparsity_opt: optimal lambda_sparsity,
  // lambda_diversity_opt: optimal lambda_diversity, 
  // cv_opt: optimal CV MSE
  
  // Initiliaze variables
  arma::uword n = x_perm.n_rows;
  arma::uword p = x_perm.n_cols;
  // Each entry of final_betas is a group of slopes, where each column
  // is the regression vector for each group
  arma::cube final_betas = zeros(p, num_groups, num_lambdas_sparsity);
  double conv_crit = 1;
  arma::uword iteration = 0;
  List ret;
  arma::vec mses_sparsity = zeros(num_lambdas_sparsity, 1);
  arma::vec mses_diversity = zeros(num_lambdas_diversity, 1);
  arma::vec lambdas_sparsity = zeros(num_lambdas_sparsity, 1);
  arma::vec lambdas_diversity = zeros(num_lambdas_diversity, 1);
  // we keep this additional copies to keep track of the last grids use,
  // for the final output
  arma::vec mses_sparsity_new = zeros(num_lambdas_sparsity, 1);
  arma::vec mses_diversity_new = zeros(num_lambdas_diversity, 1);
  arma::vec lambdas_diversity_new = zeros(num_lambdas_diversity, 1);
  
  
  double lambda_sparsity_opt = 0;
  double lambda_sparsity_opt_new = 0;
  double lambda_diversity_opt = 0;
  double lambda_diversity_opt_new = 0;
  double cv_mse_opt = 0;
  double new_mse = 0;
  // ratio of lambda_min to lambda_max   
  double eps = 0;
  
  if (n > p){
    eps = 1e-4;
  } else {
    eps = 1e-2;
  }
  
  arma::mat x_std = x_perm;
  arma::vec y_std = y_perm;
  arma::rowvec mu_x = mean(x_perm);
  arma::rowvec sd_x = stddev(x_perm, 1);
  double mu_y = mean(y_perm);
  double sd_y = stddev(y_perm, 1);
  x_std.each_row() -= mu_x;
  x_std.each_row() /= sd_x;
  y_std = y_std - mu_y;
  y_std /= sd_y;
  
  double lambda_sparsity_max = (1/ alpha) * max(abs(y_std.t() * x_std)) / n;
  lambdas_sparsity = exp(linspace(log(eps * lambda_sparsity_max), log(lambda_sparsity_max), num_lambdas_sparsity));
  bool optim_lambdas_diversity = true;
  
  while((conv_crit > tolerance) & (iteration <= 10)){
    iteration += 1;
    // Optimize over lambda_sparsity for lambda_diversity fixed
    mses_sparsity_new = CV_Ensemble_EN(x_perm, y_perm, 1, lambdas_sparsity, lambda_diversity_opt,
                                   alpha, num_groups, num_folds, tolerance,
                                   max_iter, num_threads);
    // If in first iteration, we are at the optimal
    if (iteration == 1){
      mses_sparsity = mses_sparsity_new;
      lambda_sparsity_opt = lambdas_sparsity(mses_sparsity.index_min());
      cv_mse_opt = min(mses_sparsity);
    } else {
      lambda_sparsity_opt_new = lambdas_sparsity(mses_sparsity_new.index_min());
      new_mse = min(mses_sparsity_new);
      // Check for improvement
      if (new_mse < cv_mse_opt){
        // Updating step, optimum and grid where optimum lies
        mses_sparsity = mses_sparsity_new; 
        conv_crit = (new_mse - cv_mse_opt) / new_mse;
        cv_mse_opt = new_mse;
        lambda_sparsity_opt = lambda_sparsity_opt_new;
      } else {
        // Stop everything, forget about last search
        conv_crit = 0;
        optim_lambdas_diversity = false;
      }
    }
    if(optim_lambdas_diversity){
      // Optimize over lambda_diversity for lambda_sparsity fixed
      lambdas_diversity_new = Lambdas_Diversity_Grid(x_std, y_std, lambda_sparsity_opt, num_lambdas_diversity, alpha, eps, num_groups, 
                                                 tolerance, max_iter);
      mses_diversity_new = CV_Ensemble_EN(x_perm, y_perm, 2, lambdas_diversity_new, lambda_sparsity_opt,
                                      alpha, num_groups, num_folds, tolerance,
                                      max_iter, num_threads);
      lambda_diversity_opt_new = lambdas_diversity_new[mses_diversity_new.index_min()];
      new_mse = min(mses_diversity_new);
      
      if (iteration==1){
        // in case there's no improvement
        lambdas_diversity = lambdas_diversity_new;
        mses_diversity = mses_diversity_new;
      }
      // Check for improvement
      if (new_mse < cv_mse_opt){
        // Updating step, optimum and grid where optimum lies
        lambdas_diversity = lambdas_diversity_new;
        mses_diversity = mses_diversity_new;
        conv_crit = (new_mse - cv_mse_opt) / new_mse;
        cv_mse_opt = new_mse;
        lambda_diversity_opt = lambda_diversity_opt_new;
      } else {
        // Stop everything; forget about last search
        conv_crit = 0;
      }
    } else {
      conv_crit  = 0;
    }
  }
  
  // Compute final fit using optimal tunings
  final_betas = Ensemble_EN_Grid(x_perm, y_perm, 1, lambdas_sparsity, lambda_diversity_opt,
                                 alpha, num_groups, tolerance, max_iter);
  arma::uvec index_opt = find(lambdas_sparsity == lambda_sparsity_opt);
  
  // Return elements of the list
  ret["betas"] = final_betas;
  ret["index_opt"] = index_opt[0] + 1;
  ret["lambda_sparsity_opt"] = lambda_sparsity_opt;
  ret["lambda_diversity_opt"] = lambda_diversity_opt;
  ret["lambdas_sparsity"] = lambdas_sparsity;
  ret["lambdas_diversity"] = lambdas_diversity;
  ret["cv_mse_opt"] = cv_mse_opt; 
  return(ret);
}

