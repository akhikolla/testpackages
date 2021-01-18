## \code{$summarize()} prints a summary for the fitting result under the given selector. ##
lslx$set("public",
         "summarize",
         function(selector,
                  lambda,
                  delta,
                  step,
                  standard_error = "default",
                  ridge_penalty = "default",
                  debias = "default",
                  inference = "default",
                  alpha_level = .05,
                  include_faulty = FALSE,
                  style = "default",
                  mode = "default",
                  interval = TRUE,
                  digit = 3L,
                  output = list(
                    general_information = TRUE,
                    fitting_information = FALSE,
                    saturated_model_information = FALSE,
                    baseline_model_information = FALSE,
                    numerical_condition = TRUE,
                    information_criterion = FALSE,
                    fit_index = TRUE,
                    cv_error = TRUE,
                    lr_test = TRUE,
                    rmsea_test = TRUE,
                    coefficient_test = TRUE
                  )) {
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.\n")
           }
           if (mode == "default") {
             mode <- "print"
           } else {
             if (!(mode %in% c("print", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'print' or 'return'. ")
             }
           }
           if (!(style %in% c("default", "manual", "minimal", "maximal"))) {
             stop("Argument 'style' can be only either 'default', 'manual', 'mininal', or 'maximal'.")
           }
           if (style == "default") {
             output <- list(
               general_information = TRUE,
               fitting_information = FALSE,
               saturated_model_information = FALSE,
               baseline_model_information = FALSE,
               numerical_condition = TRUE,
               information_criterion = FALSE,
               fit_index = TRUE,
               cv_error = TRUE,
               lr_test = TRUE,
               rmsea_test = TRUE,
               coefficient_test = TRUE
             )
           } else if (style == "minimal") {
             output <- 
               lapply(X = output,
                      FUN = function(output_i) {
                        output_i <- FALSE
                        return(output_i)
                      })
             output$general_information <- TRUE
             output$numerical_condition <- TRUE
           } else if (style == "maximal") {
             output <- 
               lapply(X = output,
                      FUN = function(output_i) {
                        output_i <- TRUE
                        return(output_i)
                      })
           } else if (style == "manual") {
             
           } else {}
           if (private$fitting$control$cv_fold == 1L) {
             output$cv_error <- FALSE
           }
           if (!(
             standard_error %in% c(
               "default",
               "sandwich",
               "observed_information",
               "expected_information"
             )
           )) {
             stop(
               "Argument 'standard_error' can be only either 'default', 'sandwich', 'observed_information', or 'expected_information'."
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
           ##generating output informations
           if (output$general_information) {
             general_information <-
               formatC(
                 x = c(
                   private$fitting$reduced_data$n_observation,
                   private$fitting$reduced_data$n_complete_observation,
                   ifelse(
                     private$fitting$reduced_data$n_missing_pattern == 1,
                     "none",
                     private$fitting$reduced_data$n_missing_pattern
                   ),
                   private$fitting$reduced_model$n_group,
                   private$fitting$reduced_model$n_response,
                   private$fitting$reduced_model$n_factor,
                   sum(private$fitting$reduced_model$theta_is_free),
                   sum(private$fitting$reduced_model$theta_is_pen)
                 ),
                 digits = digit,
                 format = "f"
               )
             names(general_information) <-
               c(
                 "number of observations",
                 "number of complete observations",
                 "number of missing patterns",
                 "number of groups",
                 "number of responses",
                 "number of factors",
                 "number of free coefficients",
                 "number of penalized coefficients"
               )
           } else {
             general_information <- NULL
           }
           if (output$fitting_information) {
             fitting_information <-
               formatC(
                 x = c(
                   private$fitting$control$penalty_method,
                   ifelse(
                     private$fitting$control$penalty_method %in% c("none", "forward", "backward"),
                     "none",
                     ifelse(
                       length(private$fitting$control$lambda_grid[[1]]) == 1,
                       private$fitting$control$lambda_grid[[1]],
                       paste(
                         min(private$fitting$control$lambda_grid[[1]]),
                         max(private$fitting$control$lambda_grid[[1]]),
                         sep = " - "
                       )
                     )
                   ),
                   ifelse(
                     (private$fitting$control$penalty_method %in% c("none", "forward", "backward", "lasso", "ridge")),
                     "none",
                     ifelse(
                       length(private$fitting$control$delta_grid[[1]]) == 1,
                       private$fitting$control$delta_grid[[1]],
                       paste(
                         min(private$fitting$control$delta_grid[[1]]),
                         max(private$fitting$control$delta_grid[[1]]),
                         sep = " - "
                       )
                     )
                   ),
                   ifelse(
                     (private$fitting$control$penalty_method %in% c("none", "lasso", "ridge", "elastic_net", "mcp")),
                     "none",
                     ifelse(
                       length(private$fitting$control$step_grid) == 1,
                       private$fitting$control$step_grid,
                       paste(
                         min(private$fitting$control$step_grid),
                         max(private$fitting$control$step_grid),
                         sep = " - "
                       )
                     )
                   ),
                   private$fitting$control$algorithm,
                   ifelse(
                     private$fitting$reduced_data$n_missing_pattern == 1,
                     "none",
                     sub(
                       pattern = "_",
                       replacement = " ",
                       x = private$fitting$control$missing_method
                     )
                   ),
                   private$fitting$control$tol_out
                 ),
                 digits = digit,
                 format = "f"
               )
             names(fitting_information) <-
               c(
                 "penalty method",
                 "lambda grid",
                 "delta grid",
                 "step grid",
                 "algorithm",
                 "missing method",
                 "tolerance for convergence"
               )
           } else {
             fitting_information <- NULL
           }
           
           if (output$saturated_model_information) {
             saturated_model_information <-
               formatC(
                 x = private$fitting$supplied_result$saturated_model,
                 digits = digit,
                 format = "f"
               )
             names(saturated_model_information) <-
               c("loss value",
                 "number of non-zero coefficients",
                 "degrees of freedom")
           } else {
             saturated_model_information <- NULL
           }
           
           if (output$baseline_model_information) {
             baseline_model_information <-
               formatC(
                 x = private$fitting$supplied_result$baseline_model,
                 digits = digit,
                 format = "f"
               )
             names(baseline_model_information) <-
               c("loss value",
                 "number of non-zero coefficients",
                 "degrees of freedom")
           } else {
             baseline_model_information <- NULL
           }
           
           if (output$numerical_condition) {
             numerical_condition <-
               formatC(
                 x = self$extract_numerical_condition(selector = selector,
                                                      lambda = lambda,
                                                      delta = delta,
                                                      step = step,
                                                      include_faulty = include_faulty),
                 digits = digit,
                 format = "f"
               )
             numerical_condition[["lambda_1st"]] <-
               ifelse(private$fitting$control$penalty_method %in% c("none", "forward", "backward"),
                      "none",
                      ifelse(private$fitting$control$double_regularizer,
                             paste0("c(", numerical_condition[["lambda_1st"]], ",",
                                    numerical_condition[["lambda_2nd"]], ")"),
                             numerical_condition[["lambda_1st"]]))
             numerical_condition[["delta_1st"]] <-
               ifelse(private$fitting$control$penalty_method %in% c("none", "forward", "backward", "lasso", "ridge"),
                      "none",
                      ifelse(private$fitting$control$double_regularizer,
                             paste0("c(", numerical_condition[["delta_1st"]], ",",
                                    numerical_condition[["delta_2nd"]], ")"),
                             numerical_condition[["delta_1st"]]))
             numerical_condition[["step"]] <-
               ifelse(private$fitting$control$penalty_method %in% c("none", "lasso", "ridge", "elastic_net", "mcp"),
                      "none",
                      numerical_condition[["step"]])
             numerical_condition <- numerical_condition[-c(2, 4)] 
             names(numerical_condition) <-
               c(
                 "selected lambda",
                 "selected delta",
                 "selected step",
                 "objective value",
                 "objective gradient absolute maximum",
                 "objective Hessian convexity",
                 "number of iterations",
                 "loss value",
                 "number of non-zero coefficients",
                 "degrees of freedom",
                 "robust degrees of freedom",
                 "scaling factor"
               )
           } else {
             numerical_condition <- NULL
           }
           
           if (output$information_criterion) {
             information_criterion <-
               formatC(
                 x = self$extract_information_criterion(selector = selector,
                                                        lambda = lambda,
                                                        delta = delta,
                                                        step = step,
                                                        include_faulty = include_faulty),
                 digits = digit,
                 format = "f"
               )
             names(information_criterion) <-
               c(
                 "Akaike information criterion (aic)",
                 "Akaike information criterion with penalty 3 (aic3)",
                 "consistent Akaike information criterion (caic)",
                 "Bayesian information criterion (bic)",
                 "adjusted Bayesian information criterion (abic)",
                 "Haughton Bayesian information criterion (hbic)",
                 "robust Akaike information criterion (raic)",
                 "robust Akaike information criterion with penalty 3 (raic3)",
                 "robust consistent Akaike information criterion (rcaic)",
                 "robust Bayesian information criterion (rbic)",
                 "robust adjusted Bayesian information criterion (rabic)",
                 "robust Haughton Bayesian information criterion (rhbic)"
               )
           } else {
             information_criterion <- NULL
           }
           
           if (output$fit_index) {
             fit_index <-
               formatC(
                 x = self$extract_fit_index(selector = selector,
                                             lambda = lambda,
                                             delta = delta,
                                             step = step,
                                             include_faulty = include_faulty),
                 digits = digit,
                 format = "f"
               )
             names(fit_index) <-
               c(
                 "root mean square error of approximation (rmsea)",
                 "comparative fit index (cfi)",
                 "non-normed fit index (nnfi)",
                 "standardized root mean of residual (srmr)"
               )
           } else {
             fit_index <- NULL
           }
           
           if (output$cv_error) {
             cv_error <-
               formatC(
                 x = self$extract_cv_error(selector = selector,
                                           lambda = lambda,
                                           delta = delta,
                                           step = step,
                                           include_faulty = include_faulty),
                 digits = digit,
                 format = "f"
               )
             names(cv_error) <- c("test loss value")
           } else {
             cv_error <- NULL
           }
           
           
           summary_list <-
             list(
               general_information,
               fitting_information,
               saturated_model_information,
               baseline_model_information,
               numerical_condition,
               information_criterion,
               fit_index,
               cv_error
             )
           names(summary_list) <- c(
             "General Information",
             "Fitting Information",
             "Saturated Model Information",
             "Baseline Model Information",
             "Numerical Conditions",
             "Information Criteria",
             "Fit Indices",
             "Cross-Validation Errors"
           )
           summary_list <-
             summary_list[!sapply(summary_list, is.null)]
           
           if (output$lr_test) {
             lr_test <-
               self$test_lr(selector = selector,
                            lambda = lambda,
                            delta = delta,
                            step = step,
                            include_faulty = include_faulty)
             lr_test[] <-
               data.frame(sapply(
                 X = lr_test,
                 FUN = function(lr_test_i) {
                   lr_test_i <-
                     formatC(lr_test_i, digits = digit, format = "f")
                   lr_test_i[grepl("NA", lr_test_i)] <-
                     "  -  "
                   return(lr_test_i)
                 }
               ))
             colnames(lr_test) <-
               format(c("statistic", "df", "p-value"),
                      width = 10,
                      justify = "right")
             rownames(lr_test) <-
               paste0("   ", rownames(lr_test), "  ")
           } else {
             lr_test <- NULL
           }
           
           if (output$rmsea_test) {
             rmsea_test <-
               self$test_rmsea(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 alpha_level = alpha_level,
                 include_faulty = include_faulty
               )
             rmsea_test[] <-
               data.frame(sapply(
                 X = rmsea_test,
                 FUN = function(rmsea_test_i) {
                   rmsea_test_i <-
                     formatC(rmsea_test_i, digits = digit, format = "f")
                   rmsea_test_i[grepl("NA", rmsea_test_i)] <-
                     "  -  "
                   return(rmsea_test_i)
                 }
               ))
             colnames(rmsea_test) <-
               format(colnames(rmsea_test),
                      width = 10,
                      justify = "right")
             rownames(rmsea_test) <-
               paste0("   ", rownames(rmsea_test), "  ")
           } else {
             rmsea_test <- NULL
           }
           if (output$coefficient_test) {
             coefficient_test <-
               self$test_coefficient(
                 selector = selector,
                 lambda = lambda,
                 delta = delta,
                 step = step,
                 standard_error = standard_error,
                 ridge_penalty = ridge_penalty,
                 alpha_level = alpha_level,
                 debias = debias,
                 inference = inference,
                 include_faulty = include_faulty
               )
             relation_as_groupname <-
               format(
                 private$model$specification$relation,
                 width = max(nchar(
                   private$model$specification$relation
                 )),
                 justify = "right"
               )
             block_levels <-
               c("Factor Loading",
                 "Regression",
                 "Covariance",
                 "Variance",
                 "Intercept",
                 "Threshold",
                 "Scale")
             coefficient_test[] <-
               lapply(
                 X = coefficient_test,
                 FUN = function(coefficient_test_i) {
                   coefficient_test_i <-
                     formatC(coefficient_test_i,
                             digits = digit,
                             format = "f")
                   coefficient_test_i[grepl("NA", coefficient_test_i)] <-
                     "  -  "
                   return(coefficient_test_i)
                 }
               )
             coefficient_test <-
               cbind(coefficient_test, private$model$specification[c("block",
                                                                     "group",
                                                                     "left",
                                                                     "right")])
             coefficient_test$block_label <-
               rep(NA, nrow(private$model$specification))
             coefficient_test$type <-
               format(private$model$specification$type,
                      width = 6,
                      justify = "right")
             coefficient_test$block_label[grepl("y\\|t", coefficient_test$block)] <-
               "Threshold"
             coefficient_test$block_label[grepl("(y|f)<-1", coefficient_test$block)] <-
               "Intercept"
             coefficient_test$block_label[grepl("(y|f)<-(y|f)", coefficient_test$block)] <-
               "Regression"
             coefficient_test$block_label[coefficient_test$block == "y<-f"] <-
               "Factor Loading"
             coefficient_test$block_label[grepl("(y|f)<->(y|f)", coefficient_test$block)] <-
               "Covariance"
             coefficient_test$block_label[(grepl("y<->y", coefficient_test$block) | 
                                             grepl("f<->f", coefficient_test$block)) & 
                                            (coefficient_test$left == coefficient_test$right)] <-
               "Variance"
             coefficient_test$block_label[grepl("y\\*\\*y", coefficient_test$block) ] <-
               "Scale"
             coefficient_test <-
               data.frame(coefficient_test[c(
                 "type",
                 "estimate",
                 "standard_error",
                 "z_value",
                 "p_value",
                 "lower",
                 "upper",
                 "block",
                 "block_label",
                 "group"
               )])
           } else {
             coefficient_test <- NULL
           }
           rowname_width <-
             max(nchar(unlist(lapply(
               X = summary_list, FUN = names
             )))) + 5
           value_width <-
             max(unlist(lapply(X = summary_list, FUN = nchar))) + 1
           
           ## printing summary list
           if (mode == "print") {
             for (name_i in names(summary_list)) {
               cat(name_i)
               summary_list_i <-
                 as.data.frame(summary_list[[name_i]])
               colnames(summary_list_i) <- NULL
               rownames(summary_list_i) <-
                 format(paste("  ", rownames(summary_list_i)),
                        width = rowname_width,
                        justify = "left")
               print(format(summary_list_i, width = value_width, justify = "right"))
               cat("\n")
             }
             ## printing likelihood ratio test
             if (output$lr_test) {
               cat("Likelihood Ratio Test\n")
               print(lr_test)
               cat("\n")
             }
             ## printing root mean square error of approximation test
             if (output$rmsea_test) {
               cat("Root Mean Square Error of Approximation Test\n")
               print(rmsea_test)
               cat("\n")
             }
             ## generating coefficients
             if (output$coefficient_test) {
               ## print by different groups
               if (!is.null(private$model$reference_group)) {
                 reference_group_order <-
                   which(private$model$level_group %in% private$model$reference_group)
                 group_by_order <-
                   c(reference_group_order,
                     c(1:(
                       length(private$model$level_group)
                     ))[!(c(1:(
                       length(private$model$level_group)
                     )) %in% reference_group_order)])
               } else {
                 group_by_order <- 1:length(private$model$level_group)
               }
               for (i_group in group_by_order) {
                 idc_group <-
                   coefficient_test$group == private$model$level_group[i_group]
                 coefficient_test_group <-
                   coefficient_test[idc_group, 
                                    c("type", "estimate", "standard_error",
                                      "z_value", "p_value", 
                                      "lower", "upper", "block_label"),
                                    drop = FALSE]
                 rownames(coefficient_test_group) <-
                   paste0(relation_as_groupname[idc_group], "  ")
                 colnames(coefficient_test_group) <-
                   c(
                     "type",
                     "estimate",
                     "std.error",
                     "z-value",
                     "P(>|z|)",
                     "lower",
                     "upper",
                     "block_label"
                   )
                 if (!interval) {
                   coefficient_test_group$lower <- NULL
                   coefficient_test_group$upper <- NULL
                 } 
                 if (length(private$model$level_group) != 1) {
                   cat(
                     paste0(
                       "Coefficient Test (Group = \"",
                       private$model$level_group[[i_group]],
                       "\"",
                       ", Std.Error = \"",
                       standard_error,
                       "\")\n"
                     )
                   )
                 } else {
                   cat(
                     paste0(
                       "Coefficient Test",
                       " (Std.Error = \"",
                       standard_error,
                       "\")\n"
                     )
                   )
                 }
                 ## print by block types
                 for (i_block_label in block_levels) {
                   if (sum(coefficient_test_group$block_label == i_block_label) > 0L) {
                     cat(" ", i_block_label)
                     # if 'single group' or 'reference group not specified', print nothing.
                     if ((length(group_by_order) == 1) |
                         is.null(private$model$reference_group)) {
                       cat("\n")
                     } else if (i_group == group_by_order[1]) {
                       cat(" (reference component)\n")
                     } else {
                       cat(" (increment component)\n")
                     }
                     coefficient_test_group_block <-
                       coefficient_test_group[(coefficient_test_group$block_label == i_block_label),!(colnames(coefficient_test_group) %in% c("block_label"))]
                     colnames(coefficient_test_group_block) <-
                       paste0(" ", colnames(coefficient_test_group_block))
                     print(coefficient_test_group_block)
                     cat("\n")
                   }
                 }
               }
             }
           } else if (mode == "return") {
             summary_list$'Likelihood Ratio Test' <- lr_test
             summary_list$'Root Mean Square Error of Approximation Test' <- rmsea_test
             summary_list$'Coefficient Test' <- coefficient_test
             summary_list$'Coefficient Test'$block <- NULL
             summary_list$'Coefficient Test'$block_label <- NULL
             summary_list$'Coefficient Test'$group <- NULL
             return(summary_list)
           } else {
             
           }
         })
