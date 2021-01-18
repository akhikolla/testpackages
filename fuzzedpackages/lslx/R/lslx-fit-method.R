## \code{$fit()} fits the specified model to data by minimizing a penalized ML loss function. ##
lslx$set("public",
         "fit",
         function(penalty_method = "mcp",
                  lambda_grid = "default",
                  delta_grid = "default",
                  step_grid = "default",
                  loss = "default",
                  algorithm = "default",
                  missing_method = "default",
                  start_method = "default",
                  lambda_direction = "default",
                  lambda_length = 50L,
                  delta_length = 3L,
                  threshold_value = 0.3,
                  subset = NULL,
                  cv_fold = 1L,
                  iter_out_max = 100L,
                  iter_in_max = 50L,
                  iter_other_max = 500L,
                  iter_armijo_max = 20L,
                  tol_out = 1e-3,
                  tol_in = 1e-3,
                  tol_other = 1e-7,
                  step_size = 1,
                  momentum = 0,
                  armijo = 1e-5,
                  ridge_cov = 0,
                  ridge_hessian = 1e-4,
                  ridge_weight = 1e-4,
                  warm_start = TRUE,
                  positive_variance = TRUE,
                  minimum_variance = 1e-4,
                  armijo_rule = TRUE,
                  enforce_cd = TRUE,
                  random_update = TRUE,
                  weight_matrix = NULL,
                  verbose = TRUE) {
           control <-
             list(
               penalty_method = penalty_method,
               lambda_grid = lambda_grid,
               delta_grid = delta_grid,
               step_grid = step_grid,
               loss = loss,
               algorithm = algorithm,
               missing_method = missing_method,
               start_method = start_method,
               lambda_direction = lambda_direction,
               lambda_length = lambda_length,
               delta_length = delta_length,
               threshold_value = threshold_value,
               subset = subset,
               cv_fold = cv_fold,
               iter_out_max = iter_out_max,
               iter_in_max = iter_in_max,
               iter_other_max = iter_other_max,
               iter_armijo_max = iter_armijo_max,
               tol_out = tol_out,
               tol_in = tol_in,
               tol_other = tol_other,
               step_size = step_size,
               momentum = momentum,
               armijo = armijo,
               ridge_cov = ridge_cov,
               ridge_hessian = ridge_hessian,
               warm_start = warm_start,
               positive_variance = positive_variance,
               minimum_variance = minimum_variance,
               armijo_rule = armijo_rule,
               enforce_cd = enforce_cd,
               random_update = random_update,
               weight_matrix = weight_matrix
             )
           
           private$fitting <-
             lslxFitting$new(model = private$model,
                             data = private$data,
                             control = control)
           
           if (private$fitting$control$regularizer) {
             compute_regularized_path_cpp(
               private$fitting$reduced_data,
               private$fitting$reduced_model,
               private$fitting$control,
               private$fitting$supplied_result,
               private$fitting$fitted_result
             )
           } else if (private$fitting$control$searcher) {
             compute_stepwise_path_cpp(
               private$fitting$reduced_data,
               private$fitting$reduced_model,
               private$fitting$control,
               private$fitting$supplied_result,
               private$fitting$fitted_result
             )
           } else {
             compute_none_path_cpp(
               private$fitting$reduced_data,
               private$fitting$reduced_model,
               private$fitting$control,
               private$fitting$supplied_result,
               private$fitting$fitted_result
             )
           }
           private$fitting$fitted_result$is_finite <-
             sapply(
               X = private$fitting$fitted_result$numerical_condition,
               FUN = function(numerical_condition_i) {
                 is_finite_i <- !is.na((numerical_condition_i[["objective_gradient_abs_max"]]))
                 return(is_finite_i)
               }
             )
           if (all(!private$fitting$fitted_result$is_finite)) {
             cat("ERROR: Optimization result is not finite under EVERY specified penalty level.\n",
                  "Please check model identifiability or specify better starting values.\n")
             private$fitting$fitted_result$is_convergent <-
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(numerical_condition_i) {
                   return(FALSE)
                 }
               )
             private$fitting$fitted_result$is_convex <-
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(numerical_condition_i) {
                   return(FALSE)
                 }
               )
           } else {
             private$fitting$fitted_result$is_convergent <-
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(numerical_condition_i) {
                   is_convergent_i <-
                     (numerical_condition_i[["n_iter_out"]] <= private$fitting$control$iter_out_max) &
                     (numerical_condition_i[["objective_gradient_abs_max"]] <= private$fitting$control$tol_out)
                   is_convergent_i <- 
                     ifelse(is.na(is_convergent_i), F, is_convergent_i)
                   return(is_convergent_i)
                 }
               )
             private$fitting$fitted_result$is_convex <-
               sapply(
                 X = private$fitting$fitted_result$numerical_condition,
                 FUN = function(numerical_condition_i) {
                   is_convex_i <-
                     (numerical_condition_i[["objective_hessian_convexity"]] > 0)
                   is_convex_i <- 
                     ifelse(is.na(is_convex_i), F, is_convex_i)
                   return(is_convex_i)
                 }
               )
           }
           
           if (private$fitting$control$cv_fold > 1L) {
             if (private$fitting$control$response) {
               control <- private$fitting$control
               control$start_method <- "none"       
               control$cv_fold <- 1L
               cv_error <-
                 sapply(
                   X = 1:private$fitting$control$cv_fold,
                   FUN = function(i) {
                     subset_train_i <- 
                       private$fitting$control$subset[private$fitting$control$cv_idx != i]
                     control$subset <- subset_train_i
                     fitting_train_i <-
                       lslxFitting$new(model = private$model,
                                       data = private$data,
                                       control = control)
                     fitting_train_i$supplied_result$fitted_start <- 
                       private$fitting$fitted_result$coefficient[[1]]
                     compute_regularized_path_cpp(
                       fitting_train_i$reduced_data,
                       fitting_train_i$reduced_model,
                       fitting_train_i$control,
                       fitting_train_i$supplied_result,
                       fitting_train_i$fitted_result
                     )
                     subset_test_i <-
                       private$fitting$control$subset[private$fitting$control$cv_idx == i]
                     control$subset <- subset_test_i
                     fitting_test_i <-
                       lslxFitting$new(model = private$model,
                                       data = private$data,
                                       control = control)
                     cv_error_i <-
                       sapply(
                         X = 1:length(name_grid),
                         FUN = function(j) {
                           cv_error_ij <-
                             compute_loss_value_cpp(
                               theta_value = fitting_train_i$fitted_result$coefficient[[j]],
                               reduced_data = fitting_test_i$reduced_data,
                               reduced_model = fitting_test_i$reduced_model,
                               control = fitting_test_i$control,
                               supplied_result = fitting_test_i$supplied_result)
                           return(cv_error_ij)
                         }
                       )
                   }
                 )
               private$fitting$fitted_result$cv_error <- 
                 lapply(X = apply(cv_error, 1, mean),
                        FUN = function(cv_error_i) {
                          names(cv_error_i) <- "test_loss"
                          return(cv_error_i)
                        })
             }
           }
           
           if (private$fitting$control$regularizer) {
             name_grid <-
               paste0(
                 "lambda=",
                 "c(",
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(x, "lambda_1st")
                   }
                 ),
                 ",",
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(x, "lambda_2nd")
                   }
                 ),
                 ")",
                 ",",
                 "delta=",
                 "c(",
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(x, "delta_1st")
                   }
                 ),
                 ",",
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(x, "delta_2nd")
                   }
                 ),
                 ")"
               )
           } else {
             name_grid <-
               paste0(
                 "step=",
                 sapply(
                   X = private$fitting$fitted_result$numerical_condition,
                   FUN = function(x) {
                     getElement(x, "step")
                   }
                 )
               )
           }


           names(private$fitting$fitted_result$numerical_condition) <-
             name_grid
           names(private$fitting$fitted_result$information_criterion) <-
             name_grid
           names(private$fitting$fitted_result$fit_index) <-
             name_grid
           names(private$fitting$fitted_result$cv_error) <-
             name_grid
           names(private$fitting$fitted_result$coefficient) <-
             name_grid
           names(private$fitting$fitted_result$is_finite) <- name_grid
           names(private$fitting$fitted_result$is_convergent) <- name_grid
           names(private$fitting$fitted_result$is_convex) <- name_grid
           
           if (verbose) {
             if (all(private$fitting$fitted_result$is_convergent)) {
               cat("CONGRATS: Algorithm converges under EVERY specified penalty level.\n")
               cat("  Specified Tolerance for Convergence:",
                   private$fitting$control$tol_out,
                   "\n")
               cat("  Specified Maximal Number of Iterations:",
                   private$fitting$control$iter_out_max,
                   "\n")
             } else if (all(!private$fitting$fitted_result$is_convergent)) {
               cat("WARNING: Algorithm doesn't converge under EVERY penalty level.\n")
               cat("Please try other optimization parameters or specify better starting values.\n")
             } else {
               cat("WARNING: Algorithm doesn't converge under SOME penalty level.\n")
               cat("Please try other optimization parameters or specify better starting values.\n")
             }
             if (all(!private$fitting$fitted_result$is_convex)) {
               cat("WARNING: Approximated Hessian is not convex under EVERY convexity level.\n")
               cat("Please try larger 'delta_grid' or set positive 'ridge_hessian'. \n")
             } else if (!all(private$fitting$fitted_result$is_convex)) {
               cat("WARNING: Approximated Hessian is not convex under SOME convexity level.\n")
               cat("Please try larger 'delta_grid' or set larger 'ridge_hessian'. \n")
             }
           }
         })

## \code{$fit_lasso()} fits the specified model to data by minimizing a loss function with lasso penalty (Tibshirani, 1996). ##
lslx$set("public",
         "fit_lasso",
         function(lambda_grid = "default",
                  ...) {
           self$fit(penalty_method = "lasso",
                    lambda_grid = lambda_grid,
                    ...)
         })

## \code{$fit_ridge()} fits the specified model to data by minimizing a loss function with ridge penalty (Hoerl & Kennard, 1970). ##
lslx$set("public",
         "fit_ridge",
         function(lambda_grid = "default",
                  ...) {
           self$fit(penalty_method = "ridge",
                    lambda_grid = lambda_grid,
                    ...)
         })

## \code{$fit_elastic_net()} method fits the specified model to data by minimizing a loss function with elastic net (Zou & Hastie, 2005). ##
lslx$set("public",
         "fit_elastic_net",
         function(lambda_grid = "default",
                  delta_grid = "default",
                  ...) {
           self$fit(
             penalty_method = "elastic_net",
             lambda_grid = lambda_grid,
             delta_grid = delta_grid,
             ...
           )
         })

## \code{$fit_mcp()} method fits the specified model to data by minimizing a loss function with mcp (Zhang, 2010). ##
lslx$set("public",
         "fit_mcp",
         function(lambda_grid = "default",
                  delta_grid = "default",
                  ...) {
           self$fit(
             penalty_method = "mcp",
             lambda_grid = lambda_grid,
             delta_grid = delta_grid,
             ...
           )
         })

## \code{$fit_forward()} method fits the specified model to data by minimizing a loss function with forward searching. ##
lslx$set("public",
         "fit_forward",
         function(step_grid = "default",
                  ...) {
           self$fit(
             penalty_method = "forward",
             step_grid = step_grid,
             ...
           )
         })

## \code{$fit_backward()} method fits the specified model to data by minimizing a loss function with backward searching. ##
lslx$set("public",
         "fit_backward",
         function(step_grid = "default",
                  ...) {
           self$fit(
             penalty_method = "backward",
             step_grid = step_grid,
             ...
           )
         })

## \code{$fit_none()} method fits the specified model to data by minimizing a loss function without penalty. ##
lslx$set("public",
         "fit_none",
         function(...) {
           self$fit(
             penalty_method = "none",
             ...
           )
         })