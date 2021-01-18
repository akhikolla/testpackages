## \code{$validate()} prints a summary for the validation result under the given selector.. ##
lslx$set("public",
         "validate",
         function(selector,
                  lambda,
                  delta,
                  data,
                  subset = NULL,
                  do_fit = "default",
                  standard_error = "default",
                  alpha_level = .05,
                  include_faulty = FALSE,
                  style = "default",
                  mode = "default",
                  interval = TRUE,
                  digit = 3L) {
           if (is.null(private$fitting)) {
             stop("Fitting field is not yet derived. Please use fit-related methods first.")
           }
           if (!(
             do_fit %in% c("default", "none", "level", "pattern")
           )) {
             stop(
               "Argument 'do_fit' can be only either 'default', 'none', 'level', or 'pattern'."
             )
           }
           if (do_fit == "default") {
             do_fit <- "pattern"
           } 
           if (do_fit == "none") {
             style <- "minimal"
           } 
           lslx_cv <- self$clone(deep = TRUE)
           if (!missing(data)) {
             lslx_cv$set_data(data = data)
           }
           coefficient <-
             self$extract_coefficient(selector = selector,
                                      lambda = lambda,
                                      delta = delta,
                                      include_faulty = include_faulty)
           if (do_fit %in% c("none", "pattern")) {
             type <- ifelse((coefficient != 0) & (private$model$specification$type != "fixed"),
                            "free", "fixed")
             lslx_cv$set_coefficient_start(name = names(coefficient), 
                                           start = coefficient)
             lslx_cv$set_coefficient_type(name = names(coefficient), 
                                          type = type)
             if (do_fit == "none") {
               lslx_cv$fit(penalty_method = "none",
                           algorithm = private$fitting$control$algorithm,
                           missing_method = private$fitting$control$missing_method,
                           start_method = "none",
                           subset = subset,
                           cv_fold = 1L,
                           iter_out_max = -1,
                           iter_in_max = private$fitting$control$iter_in_max,
                           iter_other_max = private$fitting$control$iter_other_max,
                           iter_armijo_max = private$fitting$control$iter_armijo_max,
                           tol_out = private$fitting$control$tol_out,
                           tol_in = private$fitting$control$tol_in,
                           tol_other = private$fitting$control$tol_other,
                           step_size = private$fitting$control$step_size,
                           armijo = private$fitting$control$armijo,
                           ridge_cov = private$fitting$control$ridge_cov,
                           ridge_hessian = private$fitting$control$ridge_hessian,
                           warm_start = private$fitting$control$warm_start,
                           positive_variance = private$fitting$control$positive_variance,
                           minimum_variance = private$fitting$control$minimum_variance,
                           enforce_cd = private$fitting$control$enforce_cd,
                           verbose = FALSE)
             } else {
               lslx_cv$fit(penalty_method = "none",
                           algorithm = private$fitting$control$algorithm,
                           missing_method = private$fitting$control$missing_method,
                           start_method = "none",
                           subset = subset,
                           cv_fold = 1L,
                           iter_out_max = private$fitting$control$iter_out_max,
                           iter_in_max = private$fitting$control$iter_in_max,
                           iter_other_max = private$fitting$control$iter_other_max,
                           iter_armijo_max = private$fitting$control$iter_armijo_max,
                           tol_out = private$fitting$control$tol_out,
                           tol_in = private$fitting$control$tol_in,
                           tol_other = private$fitting$control$tol_other,
                           step_size = private$fitting$control$step_size,
                           armijo = private$fitting$control$armijo,
                           ridge_cov = private$fitting$control$ridge_cov,
                           ridge_hessian = private$fitting$control$ridge_hessian,
                           warm_start = private$fitting$control$warm_start,
                           positive_variance = private$fitting$control$positive_variance,
                           minimum_variance = private$fitting$control$minimum_variance,
                           enforce_cd = private$fitting$control$enforce_cd,
                           verbose = FALSE)
             }
           } else if (do_fit == "level") {
               penalty_level <-
                 self$extract_penalty_level(selector = selector,
                                            lambda = lambda,
                                            delta = delta,
                                            include_faulty = include_faulty)
               lambda_ <- as.numeric(strsplit(x = penalty_level,
                                             split = "=|/")[[1]][2])
               delta_ <- as.numeric(strsplit(x = penalty_level,
                                            split = "=|/")[[1]][4])
               lslx_cv$set_coefficient_start(name = names(coefficient), 
                                             start = coefficient)
               lslx_cv$fit(penalty_method = private$fitting$control$penalty_method,
                           lambda_grid = lambda_,
                           delta_grid = delta_,
                           algorithm = private$fitting$control$algorithm,
                           missing_method = private$fitting$control$missing_method,
                           start_method = "none",
                           subset = subset,
                           cv_fold = 1L,
                           iter_out_max = private$fitting$control$iter_out_max,
                           iter_in_max = private$fitting$control$iter_in_max,
                           iter_other_max = private$fitting$control$iter_other_max,
                           iter_armijo_max = private$fitting$control$iter_armijo_max,
                           tol_out = private$fitting$control$tol_out,
                           tol_in = private$fitting$control$tol_in,
                           tol_other = private$fitting$control$tol_other,
                           step_size = private$fitting$control$step_size,
                           armijo = private$fitting$control$armijo,
                           ridge_cov = private$fitting$control$ridge_cov,
                           ridge_hessian = private$fitting$control$ridge_hessian,
                           warm_start = private$fitting$control$warm_start,
                           positive_variance = private$fitting$control$positive_variance,
                           minimum_variance = private$fitting$control$minimum_variance,
                           enforce_cd = private$fitting$control$enforce_cd,
                           verbose = FALSE)
             } else {
             }
           
           if (do_fit == "none") {
             lslx_cv$summarize(style = style,
                               mode = mode,
                               include_faulty = include_faulty,
                               interval = FALSE,
                               digit = digit)
           } else {
             lslx_cv$summarize(standard_error = standard_error,
                               alpha_level = alpha_level,
                               style = style,
                               mode = mode,
                               include_faulty = include_faulty,
                               interval = interval,
                               digit = digit)
           }
         })
