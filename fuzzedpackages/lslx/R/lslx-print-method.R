## \code{$print()} prints the \code{lslx} object. ##
lslx$set("public",
         "print",
         function() {
           cat("Brief Description of Object:\n")
           if (length(private$data$response) > 0) {
             cat("  Data field is initialized via raw data.\n")
           } else {
             cat("  Data field is initialized via moment data.\n")
           }
           if (length(private$model$name_factor) == 0) {
             if (!(c("y<-y") %in% private$model$specification$block)) {
               cat("  Model field is specified as a covariance analysis (CA) model.\n")
             } else {
               cat("  Model field is specified as a path analysis (PA) model.\n")
             }
           } else {
             if (!any(c("f<-f", "f<-y", "y<-y", "y<->f", "f<->y") %in% private$model$specification$block)) {
               cat("  Model field is specified as a factor analysis (FA) model.\n")
             } else {
               if (!any(c("f<-f", "y<->f", "f<->y") %in% private$model$specification$block)) {
                 cat("  Model field is specified as a multiple indicators and multiple causes (MIMIC) model.\n")
               } else {
                 cat("  Model field is specified as a geneal structural equation modeling (SEM) model.\n")
               }
             }
           }
           cat("    Response Variable(s):",
               private$model$name_response,
               "\n")
           if (length(private$model$name_factor) > 0) {
             cat("    Latent Factor(s):",
                 private$model$name_factor,
                 "\n") 
           }
           if (length(private$data$auxiliary) > 0) {
             cat("    Auxiliary Variable(s):",
                 colnames(private$data$auxiliary[[1]]),
                 "\n")               
           }
           if (length(private$model$level_group) > 1) {
             cat("    Group(s):",
                 private$model$level_group,
                 "\n")
             if (!is.null(private$model$reference_group)) {
               cat("    Reference Group:",
                   private$model$reference_group,
                   "\n") 
             }
           }
                  
           if (is.null(private$fitting)) {
             cat("  Fitting field is not yet derived. Please use fit-related methods.\n")
           } else {
             cat("  Fitting field is derived via maximum likelihood (ML) estimation with", 
                 private$fitting$control$penalty_method, "penalty.\n", sep = " ")
           }
           cat("\n")
           cat("Methods to Manipulate Object: \n")
           if (is.null(private$fitting)) {
             cat("  To fit the specified model to data, please use fit-related methods.\n")
             cat("    $fit() / $fit_lasso() / $fit_mcp() / $fit_none()\n")
             cat("  To modify the specified model, please use set-related methods.\n")
             cat("    $free_coefficient() / $penalize_coefficient() / $fix_coefficient()\n")
             cat("    $free_directed() / $penalize_directed() / $fix_directed()\n")
             cat("    $free_undirected() / $penalize_undirected() / $fix_undirected()\n")
             cat("    $free_block() / $penalize_block() / $fix_block()\n")
             cat("    $free_heterogeneity() / $penalize_heterogeneity() / $fix_heterogeneity()\n")
           } else {
             cat("  To summarize the fitting results, please use $summarize().\n")
             cat("  To plot the fitting results, please use plot-related methods.\n")
             cat("    $plot_numerical_condition() / $plot_information_criterion()\n")
             cat("    $plot_fit_index() / $plot_coefficients()\n")
             cat("  To obtain the test results for model or coefficients, please use test-related methods.\n")             
             cat("    $test_lr() / $test_rmsea() / $test_coefficient()\n")
             cat("  To get a deep copy of field, please use get-related methods.\n")   
             cat("    $get_model() / $get_data() / $get_fitting()\n")
             cat("  To extract specific quantity related to SEM, please use extract-related methods.\n")             
             cat("    $extract_specification() / $extract_saturated_mean() / $extract_saturated_cov()\n")
             cat("    $extract_saturated_moment_acov() \ $extract_penalty_level()\n")
             cat("    $extract_numerical_condition() / $extract_information_criterion()\n")
             cat("    $extract_fit_index() / $extract_coefficient()\n")
             cat("    $extract_implied_cov() /  $extract_implied_mean()\n")
             cat("    $extract_residual_cov() / $extract_residual_mean()\n")
             cat("    $extract_coefficient_matrix() / $extract_moment_jacobian()\n")
             cat("    $extract_expected_information() / $extract_observed_information()\n")
             cat("    $extract_score_acov() / $extract_coefficient_acov()\n")
             cat("    $extract_loss_gradient() / $extract_regularizer_gradient() / $extract_objective_gradient()\n")
           }
           cat("\n")
           cat("To learn more about 'lslx' object, please try 'help(lslx)' or see GitHub wiki (https://github.com/psyphh/lslx/wiki).")
         })
