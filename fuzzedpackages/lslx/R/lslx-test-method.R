## \code{$test_lr()} returns a \code{data.frame} of result for likelihood ratio test. ##
lslx$set("public",
         "test_lr",
         function(selector,
                  lambda,
                  delta,
                  step,
                  include_faulty = FALSE) {
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.\n")
           }
           numerical_condition <-
             self$extract_numerical_condition(selector = selector,
                                              lambda = lambda,
                                              delta = delta,
                                              step = step,
                                              include_faulty = include_faulty)
           lr_test <-
             data.frame(
               statistic = c(NA, NA),
               df = c(NA, NA),
               p_value = c(NA, NA),
               row.names = c("unadjusted", "mean-adjusted")
             )
           lr_test["unadjusted", "statistic"] <-
             numerical_condition[["loss_value"]] * private$fitting$reduced_data$n_observation
           lr_test["unadjusted", "df"] <-
             numerical_condition[["degrees_of_freedom"]]
           lr_test["unadjusted", "p_value"] <-
             1 - pchisq(lr_test["unadjusted", "statistic"],
                        lr_test["unadjusted", "df"])
           if (private$fitting$control$response) {
             scaling_factor <- numerical_condition[["scaling_factor"]]
             if (!is.na(scaling_factor)) {
               lr_test["mean-adjusted", "statistic"] <-
                 numerical_condition[["loss_value"]] * private$fitting$reduced_data$n_observation /
                 scaling_factor
               lr_test["mean-adjusted", "df"] <-
                 numerical_condition[["degrees_of_freedom"]]
               lr_test["mean-adjusted", "p_value"] <-
                 1 - pchisq(lr_test["mean-adjusted", "statistic"],
                            lr_test["mean-adjusted", "df"])
             }
           } else {
             
           }
           return(lr_test)
         })

## \code{$test_rmsea()} returns a \code{data.frame} of result for rmsea confidence intervals. ##
lslx$set("public",
         "test_rmsea",
         function(selector,
                  lambda,
                  delta,
                  step,
                  alpha_level = .05,
                  include_faulty = FALSE) {
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.\n")
           }
           numerical_condition <-
             self$extract_numerical_condition(selector = selector,
                                              lambda = lambda,
                                              delta = delta,
                                              step = step,
                                              include_faulty = include_faulty)
           fit_index <-
             self$extract_fit_index(selector = selector,
                                     lambda = lambda,
                                     delta = delta,
                                    step = step,
                                     include_faulty = include_faulty)
           lr_test <-
             self$test_lr(selector = selector,
                          lambda = lambda,
                          delta = delta,
                          step = step,
                          include_faulty = include_faulty)
           rmsea_test <-
             data.frame(
               estimate = c(NaN, NaN),
               lower = c(NaN, NaN),
               upper = c(NaN, NaN),
               row.names = c("unadjusted", "mean-adjusted")
             )
           for (row_name_i in row.names(rmsea_test)) {
             if ((row_name_i == "unadjusted") |
                 (private$fitting$control$response)) {
               lr_statistic <-
                 lr_test[row_name_i, "statistic"]
               lr_df <- lr_test[row_name_i, "df"]
               if (is.na(lr_statistic) | is.na(lr_df)) {
                 rmsea_test[row_name_i, "estimate"] <- NaN
                 rmsea_test[row_name_i, "lower"] <- NaN
                 rmsea_test[row_name_i, "upper"] <- NaN
               } else if ((lr_df == 0) &
                          (lr_statistic > sqrt(.Machine$double.eps))) {
                 rmsea_test[row_name_i, "estimate"] <- NaN
                 rmsea_test[row_name_i, "lower"] <- NaN
                 rmsea_test[row_name_i, "upper"] <- NaN
               } else if (lr_statistic < sqrt(.Machine$double.eps)) {
                 rmsea_test[row_name_i, "estimate"] <- 0
                 rmsea_test[row_name_i, "lower"] <- 0
                 rmsea_test[row_name_i, "upper"] <- 0
               } else {
                 lower_ncp <- 0
                 if (pchisq(lr_statistic,
                            lr_df, lower_ncp) < (1 - alpha_level / 2)) {
                 } else {
                   lower_ncp_1 <- lower_ncp
                   lower_ncp_2 <- 0
                   while (pchisq(lr_statistic,
                                 lr_df,
                                 lower_ncp_2) > (1 - alpha_level / 2)) {
                     lower_ncp_2 <- lower_ncp_2 + lr_df
                   }
                   lower_ncp <- (lower_ncp_1 + lower_ncp_2) / 2
                   while (abs(pchisq(lr_statistic,
                                     lr_df, lower_ncp) -
                              (1 - alpha_level / 2)) > private$fitting$control$tol_other) {
                     if (pchisq(lr_statistic,
                                lr_df,
                                lower_ncp) < (1 - alpha_level / 2)) {
                       lower_ncp_2 <- lower_ncp
                       lower_ncp <- (lower_ncp + lower_ncp_1) / 2
                     } else {
                       lower_ncp_1 <- lower_ncp
                       lower_ncp <- (lower_ncp + lower_ncp_2) / 2
                     }
                   }
                 }
                 upper_ncp <- 0
                 if (pchisq(lr_statistic,
                            lr_df, upper_ncp) < (alpha_level / 2)) {
                 } else {
                   upper_ncp_1 <- upper_ncp
                   upper_ncp_2 <- 0
                   while (pchisq(lr_statistic,
                                 lr_df,
                                 upper_ncp_2) > (alpha_level / 2)) {
                     upper_ncp_2 <- upper_ncp_2 + lr_df
                   }
                   upper_ncp <- (upper_ncp_1 + upper_ncp_2) / 2
                   while (abs(pchisq(lr_statistic,
                                     lr_df, upper_ncp) -
                              (alpha_level / 2)) > private$fitting$control$tol_other) {
                     if (pchisq(lr_statistic,
                                lr_df,
                                upper_ncp) < (alpha_level / 2)) {
                       upper_ncp_2 <- upper_ncp
                       upper_ncp <- (upper_ncp + upper_ncp_1) / 2
                     } else {
                       upper_ncp_1 <- upper_ncp
                       upper_ncp <- (upper_ncp + upper_ncp_2) / 2
                     }
                   }
                 }
                 if (row_name_i == "unadjusted") {
                   scaling_factor <- 1
                 } else {
                   scaling_factor <- numerical_condition[["scaling_factor"]]
                 }
                 if (!is.na(scaling_factor)) {
                   rmsea_test[row_name_i, "estimate"] <-
                     sqrt(max(
                       0,
                       scaling_factor * private$fitting$reduced_model$n_group * (lr_statistic - lr_df) /
                         (private$fitting$reduced_data$n_observation * lr_df)
                     ))
                   rmsea_test[row_name_i, "lower"]  <-
                     sqrt(
                       max(
                         0,
                         scaling_factor * private$fitting$reduced_model$n_group * lower_ncp /
                           (private$fitting$reduced_data$n_observation * lr_df)
                       )
                     )
                   rmsea_test[row_name_i, "upper"]  <-
                     sqrt(
                       max(
                         0,
                         scaling_factor * private$fitting$reduced_model$n_group * upper_ncp /
                           (private$fitting$reduced_data$n_observation * lr_df)
                       )
                     )
                 }
               }
             }
           }
           return(rmsea_test)
         })

## \code{$test_coefficient()} returns a \code{data.frame} of result for coefficient significance and confidence interval. ##
lslx$set("public",
         "test_coefficient",
         function(selector,
                  lambda,
                  delta,
                  step,
                  standard_error = "default",
                  ridge_penalty = "default",
                  debias = "default",
                  inference = "default",
                  alpha_level = .05,
                  include_faulty = FALSE) {
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.\n")
           }
           if (!(
             standard_error %in% c("default", "sandwich", "observed_information", "expected_information")
           )) {
             stop(
               "Argument 'standard_error' can be only either 'default', 'sandwich', 'observed_information', or 'expected_information'."
             )
           }
           if (!(
             inference %in% c("default", "naive", "polyhedral", "scheffe")
           )) {
             stop(
               "Argument 'inference' can be only either 'default', 'naive', 'polyhedral', or 'scheffe'."
             )
           }
           if (!(
             debias %in% c("default", "none", "one_step")
           )) {
             stop(
               "Argument 'debias' can be only either 'default', 'none', or 'one_step'."
             )
           }
           if (standard_error == "default") {
             if (private$fitting$control$response) {
               standard_error <- "sandwich"
             } else {
               standard_error <- "observed_information"
             }
           }
           if (ridge_penalty == "default") {
             if (private$fitting$control$penalty_method %in% c("ridge", "elastic_net")) {
               ridge_penalty <- TRUE
             } else {
               ridge_penalty <- FALSE
             }
           }
           if (inference == "default") {
             inference <- "naive"
             if (debias == "default") {
               debias <- "none"
             }
           } else if (inference == "polyhedral") {
             if (debias == "default") {
               debias <- "one_step"
             }
             if (debias == "none") {
               stop("Argument 'debias' cannot be 'none' under 'inference' == 'polyhedral'.")
             }
           } else if (inference == "scheffe") {
             if (debias == "default") {
               debias <- "none"
             }
           } else {}
           if (private$fitting$control$penalty_method == "none") {
             if (inference %in% c("polyhedral", "scheffe")) {
               stop(
                 "When 'penalty_method' is 'none', 'inference' cannot be 'polyhedral' or 'scheffe'."
               )
             }
           }
           if (private$fitting$control$penalty_method %in% c("forward", "backward", "ridge")) {
             if (inference %in% c("polyhedral")) {
               stop(
                 "When 'penalty_method' is 'forward' or 'backward', 'inference' cannot be 'polyhedral'."
               )
             }
           }
           
           coefficient <-
             self$extract_coefficient(selector = selector,
                                      lambda = lambda,
                                      delta = delta,
                                      step = step,
                                      include_faulty = include_faulty)
           coefficient_acov <-
             self$extract_coefficient_acov(
               selector = selector,
               lambda = lambda,
               delta = delta,
               step = step,
               standard_error = standard_error,
               ridge_penalty = ridge_penalty,
               include_faulty = include_faulty
             )
           if (debias == "none") {
             coefficient_test <-
               data.frame(estimate = coefficient,
                          standard_error = sqrt(diag(coefficient_acov)))             
           } else {
             debiased_coefficient <-
               self$extract_debiased_coefficient(selector = selector,
                                                 lambda = lambda,
                                                 delta = delta,
                                                 step = step,
                                                 include_faulty = include_faulty)
             coefficient_test <-
               data.frame(estimate = debiased_coefficient,
                          standard_error = sqrt(diag(coefficient_acov)))
           }
           coefficient_test$z_value <-
             coefficient_test$estimate / coefficient_test$standard_error
           attr(coefficient_test, "standard_error") <-
             standard_error
           if (inference == "naive") {
             coefficient_test$p_value <-
               2 * pnorm(-abs(coefficient_test$z_value))
             coefficient_test$lower <-
               coefficient_test$estimate + qnorm(alpha_level / 2) * coefficient_test$standard_error
             coefficient_test$upper <-
               coefficient_test$estimate + qnorm(1 - alpha_level / 2) * coefficient_test$standard_error
           } else {
             is_active <-
               private$fitting$reduced_model$theta_is_free |
               (private$fitting$reduced_model$theta_is_pen &
                  coefficient != 0)
             is_pen <- private$fitting$reduced_model$theta_is_pen
             is_selected <- is_pen & (coefficient != 0)
             if (!any(is_selected)) {
               stop(
                 "No non-zero parameters are selected and hence post-selection inference cannot be implemented."
               )
             }
             
             if (inference == "polyhedral") {
               a_ph <- - diag(sign(coefficient))
               b_ph <-
                 matrix((sign(coefficient) * (coefficient - debiased_coefficient)))
               tnorm_inference <-
                 lapply(
                   X = 1:length(debiased_coefficient),
                   FUN = function(i) {
                     tnorm_quantity <- 
                       compute_tnorm_quantity(i, a_ph, b_ph, 
                                              debiased_coefficient, coefficient_acov,
                                              is_pen, is_active, is_selected) 
                     tnorm_p_value <-  
                       compute_tnorm_p_value(theta = tnorm_quantity$theta, 
                                             mu = 0, sigma = tnorm_quantity$sigma, 
                                             left = tnorm_quantity$left, 
                                             right = tnorm_quantity$right)
                     tnorm_interval <- 
                       compute_tnorm_interval(theta = tnorm_quantity$theta, 
                                              sigma = tnorm_quantity$sigma, 
                                              left = tnorm_quantity$left, 
                                              right = tnorm_quantity$right, 
                                              alpha_level = alpha_level, 
                                              grid_range = c(-100, 100), 
                                              grid_length = 100, 
                                              depth_max = 3) 
                     return(c(p_value = tnorm_p_value,
                              lower = tnorm_interval[["lower"]],
                              upper = tnorm_interval[["upper"]]))
                   }
                 )
               coefficient_test$p_value <- 
                 sapply(X = tnorm_inference, 
                        FUN = function(tnorm_inference_i) {
                          getElement(tnorm_inference_i, "p_value")
                        })
               coefficient_test$lower <-
                 sapply(X = tnorm_inference, 
                        FUN = function(tnorm_inference_i) {
                          getElement(tnorm_inference_i, "lower")
                        })
               coefficient_test$upper <- 
                 sapply(X = tnorm_inference, 
                        FUN = function(tnorm_inference_i) {
                          getElement(tnorm_inference_i, "upper")
                        })
             }  else if (inference == "scheffe") {
               df_scheffe <- sum(private$fitting$reduced_model$theta_is_pen)
               c_scheffe <- sqrt(qchisq(1 - alpha_level, df = df_scheffe))
               coefficient_test$p_value <-
                 ifelse(private$fitting$reduced_model$theta_is_pen,
                        1 - pchisq((coefficient_test$z_value)^2, df = df_scheffe),
                        2 * pnorm(-abs(coefficient_test$z_value)))
               coefficient_test$lower <-
                 ifelse(private$fitting$reduced_model$theta_is_pen,
                        coefficient_test$estimate - c_scheffe * coefficient_test$standard_error,
                        coefficient_test$estimate + qnorm(alpha_level / 2) * coefficient_test$standard_error)
                 
               coefficient_test$upper <-
                 ifelse(private$fitting$reduced_model$theta_is_pen,
                        coefficient_test$estimate + c_scheffe * coefficient_test$standard_error,
                        coefficient_test$estimate + qnorm(1 - alpha_level / 2) * coefficient_test$standard_error)
             }
           }
           return(coefficient_test)
         })


