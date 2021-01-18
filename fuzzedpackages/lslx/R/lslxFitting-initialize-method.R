## \code{$new()} initializes a new \code{lslxFitting} object. ##
lslxFitting$set("public",
                "initialize",
                function(model,
                         data,
                         control) {
                  private$initialize_control(model = model,
                                             data = data,
                                             control = control)
                  private$initialize_reduced_model(model = model)
                  private$initialize_reduced_data(data = data)
                  private$initialize_weight_matrix()
                  self$supplied_result <- list()
                  private$compute_baseline_model()
                  private$compute_saturated_model()
                  private$initialize_supplied_result()
                  private$initialize_grid()
                  private$initialize_fitted_result()
                })


## \code{$initialize_control()} initializes control options. ##
lslxFitting$set("private",
                "initialize_control",
                function(model,
                         data,
                         control) {
                  self$control <- control
                  if (!(
                    self$control$penalty_method %in% c(
                      "none",
                      "lasso",
                      "ridge",
                      "elastic_net",
                      "mcp",
                      "forward",
                      "backward"
                    )
                  )) {
                    stop(
                      "Argument 'penalty_method' can be only either 'none', 'lasso', 'ridge', 'elastic_net', 'mcp', 'forward', or 'backward'."
                    )
                  }
                  if (!is.numeric(self$control$lambda_grid)) {
                    if (!is.character(self$control$lambda_grid)) {
                      stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$lambda_grid) &
                               (length(self$control$lambda_grid) != 1)) {
                      stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$lambda_grid) &
                               (length(self$control$lambda_grid) == 1)) {
                      if (self$control$lambda_grid != "default") {
                        stop("Argument 'lambda_grid' can be only a numeric vector or set as 'default'.")
                      }
                    }
                  }
                  if (!is.numeric(self$control$delta_grid)) {
                    if (!is.character(self$control$delta_grid)) {
                      stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$delta_grid) &
                               (length(self$control$delta_grid) != 1)) {
                      stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$delta_grid) &
                               (length(self$control$delta_grid) == 1)) {
                      if (self$control$delta_grid != "default") {
                        stop("Argument 'delta_grid' can be only a numeric vector or set as 'default'.")
                      }
                    }
                  }
                  if (!is.numeric(self$control$step_grid)) {
                    if (!is.character(self$control$step_grid)) {
                      stop("Argument 'step_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$step_grid) &
                               (length(self$control$step_grid) != 1)) {
                      stop("Argument 'step_grid' can be only a numeric vector or set as 'default'.")
                    } else if (is.character(self$control$step_grid) &
                               (length(self$control$step_grid) == 1)) {
                      if (self$control$step_grid != "default") {
                        stop("Argument 'step_grid' can be only a numeric vector or set as 'default'.")
                      }
                    }
                  }
                  
                  if (!(self$control$loss %in% c("default", "ml", "uls", "dwls", "wls"))) {
                    stop("Argument 'default' can be only 'default', 'ml', 'uls', 'dwls', or 'wls'.")
                  }
                  if (!(self$control$algorithm %in% c("default", "dynamic", "gd", "bfgs", "fisher"))) {
                    stop("Argument 'algorithm' can be only 'default', 'dynamic', 'gd', 'bfgs', or 'fisher'.")
                  }
                  if (!(self$control$missing_method %in% c("default", "two_stage", "listwise_deletion"))) {
                    stop(
                      "Argument 'missing_method' can be only 'default', 'two_stage', or 'listwise_deletion'."
                    )
                  }
                  if (!(self$control$start_method %in% c("default", "none", "mh", "heuristic"))) {
                    stop("Argument 'start_method' can be only 'default', 'none', 'mh', or 'heuristic'.")
                  }
                  if (!(self$control$lambda_direction %in% c("default", "manual", "decrease", "increase"))) {
                    stop(
                      "Argument 'lambda_direction' can be only 'default', 'manual', 'decrease', or 'increase'."
                    )
                  }
                  if (!is.null(self$control$subset)) {
                    if (!(is.integer(self$control$subset) |
                          is.logical(self$control$subset))) {
                      stop("Argument 'subset' must be a integer or logical vector.")
                    }
                  }
                  if (!(is.numeric(self$control$lambda_length) &
                        (self$control$lambda_length > 0))) {
                    stop("Argument 'lambda_length' must be a positive integer.")
                  }
                  if (!(is.numeric(self$control$delta_length) &
                        (self$control$delta_length > 0))) {
                    stop("Argument 'delta_length' must be a positive integer.")
                  }
                  if (!(is.numeric(self$control$threshold_value) &
                        (self$control$threshold_value > 0))) {
                    stop("Argument 'threshold_value' must be a positive value.")
                  }
                  if (!(is.numeric(self$control$cv_fold) &
                        (self$control$cv_fold > 0))) {
                    stop("Argument 'cv_fold' must be a positive integer.")
                  }
                  if (!(is.numeric(self$control$iter_out_max) &
                        (length(self$control$iter_out_max) == 1))) {
                    stop("Argument 'iter_out_max' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$iter_in_max) &
                        (length(self$control$iter_in_max) == 1))) {
                    stop("Argument 'iter_in_max' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$iter_armijo_max) &
                        (length(self$control$iter_armijo_max) == 1))) {
                    stop("Argument 'iter_armijo_max' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$tol_out) &
                        (length(self$control$tol_out) == 1))) {
                    stop("Argument 'tol_out' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$tol_in) &
                        (length(self$control$tol_in) == 1))) {
                    stop("Argument 'tol_in' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$step_size) &
                        (length(self$control$step_size) == 1))) {
                    stop("Argument 'step_size' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$armijo) &
                        (length(self$control$armijo) == 1))) {
                    stop("Argument 'armijo' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$ridge_cov) &
                        (length(self$control$ridge_cov) == 1))) {
                    stop("Argument 'ridge_cov' must be a numeric vector with length one.")
                  }
                  if (!(is.numeric(self$control$ridge_hessian) &
                        (length(self$control$ridge_hessian) == 1))) {
                    stop("Argument 'ridge_hessian' must be a numeric vector with length one.")
                  }
                  if (!(is.logical(self$control$warm_start) &
                        (length(self$control$warm_start) == 1))) {
                    stop("Argument 'warm_start' must be a logical vector with length one.")
                  }
                  if (!(is.logical(self$control$positive_variance) &
                        (length(self$control$positive_variance) == 1))) {
                    stop("Argument 'positive_variance' must be a logical vector with length one.")
                  }
                  if (!(is.numeric(self$control$minimum_variance) &
                        (length(self$control$minimum_variance) == 1))) {
                    stop("Argument 'minimum_variance' must be a numeric vector with length one.")
                  }
                  if (!(is.logical(self$control$enforce_cd) &
                        (length(self$control$enforce_cd) == 1))) {
                    stop("Argument 'enforce_cd' must be a logical vector with length one.")
                  }
                  if (!(is.logical(self$control$random_update) &
                        (length(self$control$random_update) == 1))) {
                    stop("Argument 'random_update' must be a logical vector with length one.")
                  }
                  if (!is.null(self$control$weight_matrix)) {
                    if (!is.list(self$control$weight_matrix) &
                        !is.matrix(self$control$weight_matrix)) {
                      stop(
                        "Argument 'weight_matrix' must be a 'matrix' (for single group analysis)",
                        " or a 'list' of 'matrix' (for multiple group analysis)."
                      )
                    }
                    if (is.matrix(self$control$weight_matrix)) {
                      self$control$weight_matrix <- list(self$control$weight_matrix)
                    }
                    if (length(self$control$weight_matrix) != length(model$level_group)) {
                      stop(
                        "The length of argument 'weight_matrix' doesn't match the number of groups.",
                        "\n  The length of 'weight_matrix' is ",
                        length(self$control$weight_matrix),
                        ".",
                        "\n  The number of groups is ",
                        length(model$level_group),
                        "."
                      )
                    }
                    if (length(dim(self$control$weight_matrix[[1]])) != 2 |
                        !is.numeric(self$control$weight_matrix[[1]])) {
                      stop("Some element in 'weight_matrix' is not a matrix.")
                    }
                  }
                  if (length(model$name_factor) == 0) {
                    if (!(c("y<-y") %in% model$specification$block)) {
                      self$control$model_class <- "ca"
                    } else {
                      self$control$model_class <- "pa"
                    }
                  } else {
                    if (!any(c("f<-f", "f<-y", "y<-y", "y<->f", "f<->y") %in% model$specification$block)) {
                      if ("y<-f" %in% model$specification$block) {
                        if (all(model$specification$type[model$specification$block == "y<-f"] == "pen")) {
                          self$control$model_class <- "efa"
                        } else {
                          self$control$model_class <- "fa"
                        }
                      } else {
                        self$control$model_class <- "sem"
                      }
                    } else {
                      if (!any(c("f<-f", "y<->f", "f<->y") %in% model$specification$block)) {
                        self$control$model_class <- "mimic"
                      } else {
                        self$control$model_class <- "sem"
                      }
                    }
                  }
                  if (2 %in% model$specification$set) {
                    self$control$double_regularizer <- TRUE
                  } else {
                    self$control$double_regularizer <- FALSE
                  }
                  if (length(data$response) > 0) {
                    self$control$response <- TRUE
                  } else {
                    self$control$response <- FALSE
                  }
                  if (length(data$auxiliary) > 0) {
                    self$control$auxiliary <- TRUE
                  } else {
                    self$control$auxiliary <- FALSE
                  }
                  if (length(model$ordered_variable) == 0) {
                    self$control$continuous <- TRUE
                  } else {
                    self$control$continuous <- FALSE
                  }
                  if (self$control$penalty_method %in% c("lasso", "ridge", "elastic_net", "mcp")) {
                    self$control$regularizer <- TRUE
                  } else {
                    self$control$regularizer <- FALSE
                  }
                  if (self$control$penalty_method %in% c("forward", "backward")) {
                    self$control$searcher <- TRUE
                  } else {
                    self$control$searcher <- FALSE
                  }
                  if (self$control$penalty_method == "none" &
                      any(model$specification$type == "pen")) {
                    stop(
                      "When the specified model includes penalized coefficients, 'penalty_method' cannot be 'none'."
                    )
                  }
                  if (self$control$regularizer) {
                    self$control$regularizer_type <- self$control$penalty_method
                    if (self$control$lambda_grid[[1]] != "default") {
                      if (any(self$control$lambda_grid < 0)) {
                        stop(
                          "When 'penalty_method' is set as 'lasso' or 'mcp', any element in 'lambda_grid' must be non-negative."
                        )
                      }
                    }
                    if (self$control$delta_grid[[1]] != "default") {
                      if (self$control$penalty_method %in% c("lasso", "ridge")) {
                        stop(
                          "When 'penalty_method' is set as 'lasso' or 'ridge', 'delta_grid' must be set as 'default'."
                        )
                      } else if (self$control$penalty_method == "elastic_net") {
                        if (any(self$control$delta_grid < 0) |
                            any(self$control$delta_grid > 1)) {
                          stop(
                            "When 'penalty_method' is set as 'elastic_net', any element in 'delta_grid' must be in [0, 1]."
                          )
                        }
                      } else if (self$control$penalty_method == "mcp") {
                        if (any(self$control$delta_grid <= 0)) {
                          stop(
                            "When 'penalty_method' is set as 'mcp', any element in 'delta_grid' must be positive."
                          )
                        }
                      } else {
                      }
                    }
                  } else {
                    self$control$regularizer_type <- "none"
                    if (self$control$lambda_grid[[1]] != "default") {
                      stop(
                        "When 'penalty_method' is set as 'none', 'forward', or 'backward', 'lambda_grid' must be set as 'default'."
                      )
                    }
                    if (self$control$delta_grid[[1]] != "default") {
                      stop(
                        "When 'penalty_method' is set as 'none', 'forward', or 'backward', 'delta_grid' must be set as 'default'."
                      )
                    }
                    if (self$control$lambda_direction != "default") {
                      stop(
                        "When 'penalty_method' is set as 'none', 'forward', or 'backward', 'lambda_direction' must be set as 'default'."
                      )
                    }
                  }
                  
                  if (self$control$searcher) {
                    self$control$searcher_type <- self$control$penalty_method
                  } else {
                    self$control$searcher_type <- "none"
                    if (self$control$step_grid[[1]] != "default") {
                      stop(
                        "When 'penalty_method' is set as 'lasso', 'ridge', 'elastic_net', or mcp', 'step_grid' must be set as 'default'."
                      )
                    }
                  }
                  
                  if (self$control$loss == "default") {
                    if (self$control$continuous) {
                      self$control$loss <- "ml"
                    } else {
                      self$control$loss <- "dwls"
                    }
                  } else {
                    if (!self$control$continuous) {
                      if (self$control$loss == "ml") {
                        stop("When ordered variable is included in the model, 'loss' cannot be 'ml'.")
                      } 
                    } 
                  }
                  
                  if (!is.null(self$control$weight_matrix)) {
                    if (self$control$loss != "wls") {
                      stop("When 'weight_matrix' is specified, 'loss' must be set as 'wls'.")
                    }
                  }
                  if (self$control$algorithm == "default") {
                    self$control$algorithm <- "dynamic"
                  }
                  if (self$control$missing_method == "default") {
                    if (self$control$continuous) {
                      if (self$control$response) {
                        self$control$missing_method <- "two_stage"
                      } else {
                        self$control$missing_method <- "listwise_deletion"
                      }
                    } else {
                      self$control$missing_method <- "listwise_deletion"
                    }
                  }
                  if (self$control$start_method == "default") {
                    self$control$start_method <- "mh"
                  }
                  if (self$control$regularizer &
                      (!self$control$enforce_cd)) {
                    stop("When any regularizer is specified, 'enforce_cd' must be set as TRUE.")
                  }
                  if (!is.null(self$control$subset)) {
                    if (self$control$response) {
                      if (is.logical(self$control$subset)) {
                        self$control$subset <- which(self$control$subset)
                      }
                    } else {
                      stop("When only moment data is available, 'subset' must be set as NULL.")
                    }
                  } else {
                    if (self$control$response) {
                      self$control$subset <-
                        1:sum(sapply(data$response, FUN = nrow))
                    } else {
                      self$control$subset <-
                        1:sum(unlist(data$sample_size))
                    }
                  }
                  if (!is.integer(self$control$cv_fold)) {
                    self$control$cv_fold <- as.integer(self$control$cv_fold)
                  }
                  if (self$control$cv_fold > 1L) {
                    if (self$control$response) {
                      cv_idx <-
                        sample(
                          x = 1:self$control$cv_fold,
                          size = length(self$control$subset),
                          replace = TRUE
                        )
                      self$control$cv_idx <- cv_idx
                    }
                  }
                })

## \code{$initialize_reduced_model()} initializes a reduced model. ##
lslxFitting$set("private",
                "initialize_reduced_model",
                function(model) {
                  self$reduced_model <-
                    list(
                      n_response = length(model$name_response),
                      n_factor =  length(model$name_factor),
                      n_eta = length(model$name_eta),
                      n_moment = ifelse(self$control$continuous,
                                        length(model$name_response) *
                                          (length(model$name_response) + 3) / 2,
                                        length(model$numeric_variable) * (length(model$numeric_variable)  + 1) / 2 +
                                          length(model$ordered_variable) * (length(model$ordered_variable) - 1) / 2 +
                                          length(model$numeric_variable) * length(model$ordered_variable)+
                                          length(model$numeric_variable) + sum(model$nlevel_ordered - 1)),
                      n_moment_1 = ifelse(self$control$continuous,
                                        length(model$name_response),
                                          length(model$numeric_variable) + sum(model$nlevel_ordered - 1)),
                      n_moment_2 = ifelse(self$control$continuous,
                                        length(model$name_response) *
                                          (length(model$name_response) + 1) / 2,
                                        length(model$numeric_variable) * (length(model$numeric_variable)  + 1) / 2 +
                                          length(model$ordered_variable) * (length(model$ordered_variable) - 1) / 2 +
                                          length(model$numeric_variable) * length(model$ordered_variable)),
                      n_group = length(model$level_group),
                      n_numeric = length(model$numeric_variable),
                      n_ordered = length(model$ordered_variable),
                      n_threshold = sum(model$nlevel_ordered - 1),
                      n_theta = nrow(model$specification),
                      n_theta_is_free = NA_integer_,
                      n_theta_is_pen = NA_integer_,
                      idx_numeric = match(model$numeric_variable, model$name_response),
                      idx_ordered = match(model$ordered_variable, model$name_response),
                      idx_reference = ifelse(is.null(model$reference_group), 
                                             0, 
                                             ifelse(is.na(model$reference_group), 
                                                    0, 
                                                    match(model$reference_group, model$level_group))),
                      eta_is_exogenous = model$name_eta %in% model$name_exogenous,
                      eta_is_endogenous = model$name_eta %in% model$name_endogenous,
                      theta_name = rownames(model$specification),
                      theta_matrix_idx =
                        ifelse(
                          model$specification$matrix == "gamma",
                          0L,
                          ifelse(
                            model$specification$matrix == "alpha",
                            1L,
                            ifelse(model$specification$matrix == "beta",
                                   2L,
                                   ifelse(model$specification$matrix == "phi",
                                          3L,
                                          4L))
                          )
                        ),
                      theta_left_idx = match(model$specification$left, model$name_eta),
                      theta_right_idx = ifelse(
                        model$specification$matrix == "gamma",
                        match(model$specification$right, model$name_threshold),
                        ifelse(
                          model$specification$matrix == "alpha",
                          1L,
                          match(model$specification$right, model$name_eta)
                        )
                      ),
                      theta_flat_idx = NA_integer_,
                      theta_tflat_idx = NA_integer_,
                      theta_group_idx = ifelse(
                        rep(
                          is.null(model$reference_group),
                          length(model$specification$group)
                        ),
                        match(model$specification$group,
                              model$level_group),
                        ifelse(
                          model$specification$group ==
                            model$reference_group,
                          0L,
                          match(model$specification$group,
                                model$level_group)
                        )
                      ),
                      theta_is_free = (model$specification$type == "free"),
                      theta_is_pen = (model$specification$type == "pen"),
                      theta_is_diag =
                        ifelse(
                          model$specification$matrix == "phi" &
                            (model$specification$left ==
                               model$specification$right),
                          TRUE,
                          FALSE
                        ),
                      theta_start = model$specification$start,
                      theta_penalty = model$specification$penalty,
                      theta_set = model$specification$set,
                      theta_weight = model$specification$weight)
                  self$reduced_model$theta_penalty <-
                    ifelse(self$reduced_model$theta_penalty == "default",
                           self$control$penalty_method,
                           self$reduced_model$theta_penalty)
                  if (self$control$penalty_method %in% c("forward", "backward")) {
                    if (any(self$reduced_model$theta_penalty %in% c("lasso", "ridge", "mcp", "elastic_net"))) {
                      stop("When 'penalty_method' is set as 'forward' or 'backward', coefficients cannot be penalized by 'lasso', 'ridge', 'mcp', or 'elastic_net'.")
                    }
                  }

                  if (!self$control$continuous) {
                    self$reduced_model$threshold_flat_idx <- 
                      matrix(0, length(model$name_response), max(model$nlevel_ordered - 1))
                    self$reduced_model$idx_gamma <- integer()
                    self$reduced_model$idx_mu <- integer()
                    self$reduced_model$idx_sigma <- integer()
                    count_threshold <- 1L
                    for (i in 1:length(model$name_response)) {
                      if (model$name_response[i] %in% model$ordered_variable) {
                        nlevel_ordered_i <- 
                          model$nlevel_ordered[which(model$name_response[i] == model$ordered_variable)]
                        for (j in 1:(nlevel_ordered_i - 1L)) {
                          self$reduced_model$threshold_flat_idx[i, j] <- 
                            count_threshold
                          count_threshold <- count_threshold + 1L
                        }
                      } 
                    }
                    
                    count_threshold <- 1L
                    for (response_i in model$name_response) {
                      if (response_i %in% model$ordered_variable) {
                        nlevel_ordered_i <- 
                          model$nlevel_ordered[which(response_i == model$ordered_variable)]
                        for (j in 1:(nlevel_ordered_i - 1L)) {
                          self$reduced_model$idx_gamma <- 
                            c(self$reduced_model$idx_gamma, count_threshold)
                          count_threshold <- count_threshold + 1L
                        }
                      } else {
                        self$reduced_model$idx_gamma <- 
                          c(self$reduced_model$idx_gamma, 0L)
                      }
                    }
                    
                    for (response_i in model$name_response) {
                      if (response_i %in% model$ordered_variable) {
                        nlevel_ordered_i <- 
                          model$nlevel_ordered[which(response_i == model$ordered_variable)]
                        for (j in 1:(nlevel_ordered_i - 1L)) {
                          self$reduced_model$idx_mu <- 
                            c(self$reduced_model$idx_mu, 
                              which(response_i == model$name_response))
                        }
                      } else {
                        self$reduced_model$idx_mu <- 
                          c(self$reduced_model$idx_mu, 
                            which(response_i == model$name_response))
                      }
                    }
                    
                    for (i in 1:self$reduced_model$n_response) {
                      if (i %in% self$reduced_model$idx_numeric) {
                        self$reduced_model$idx_sigma <- 
                          c(self$reduced_model$idx_sigma, 
                            as.integer(self$reduced_model$n_response * (i-1) + 
                                         i - i * (i - 1L) / 2L))
                      }
                    }
                    if (self$reduced_model$n_response > 1) {
                      for (i in 1:(self$reduced_model$n_response - 1)) {
                        for (j in (i+1):self$reduced_model$n_response) {
                          self$reduced_model$idx_sigma <- 
                            c(self$reduced_model$idx_sigma, 
                              as.integer(self$reduced_model$n_response * (i-1) + 
                                           j - i * (i - 1L) / 2L) )
                        }
                      }  
                    }
                    self$reduced_model$idx_diag <- 
                      diag(matrix(1:(length(model$name_eta) * length(model$name_eta)), 
                                  length(model$name_eta), length(model$name_eta)))[1:length(model$name_response)]
                    self$reduced_model$idx_nondiag <- 
                      setdiff(matrix(1:(length(model$name_eta) * length(model$name_eta)), 
                                     length(model$name_eta), length(model$name_eta)),
                              self$reduced_model$idx_diag)
                    self$reduced_model$idx_diag_psi <- 
                      1L + ((1:self$reduced_model$n_response) - 1L) * (self$reduced_model$n_response + 1L)
                  }
                  
                  self$reduced_model$theta_flat_idx <-
                    ifelse(self$reduced_model$theta_matrix_idx == 0,
                           self$reduced_model$threshold_flat_idx[
                             cbind(ifelse(self$reduced_model$theta_left_idx <= 
                                            length(model$name_response),
                                          self$reduced_model$theta_left_idx,
                                          NA),
                                   ifelse(self$reduced_model$theta_right_idx <=
                                            length(model$name_threshold),
                                          self$reduced_model$theta_right_idx,
                                          NA))],
                           ifelse(
                             self$reduced_model$theta_matrix_idx == 1,
                             self$reduced_model$theta_left_idx,
                             ifelse(
                               self$reduced_model$theta_matrix_idx == 2,
                               self$reduced_model$n_eta *
                                 (self$reduced_model$theta_right_idx - 1L) +
                                 self$reduced_model$theta_left_idx,
                               ifelse(self$reduced_model$theta_matrix_idx == 3,
                                      as.integer(
                                        self$reduced_model$n_eta *
                                          (self$reduced_model$theta_right_idx - 1L) +
                                          self$reduced_model$theta_left_idx
                                      ),
                                      as.integer(self$reduced_model$theta_left_idx))
                             )
                           ))
                  self$reduced_model$theta_tflat_idx <-
                    ifelse(self$reduced_model$theta_matrix_idx == 0,
                           self$reduced_model$threshold_flat_idx[
                             cbind(ifelse(self$reduced_model$theta_left_idx <= 
                                            length(model$name_response),
                                          self$reduced_model$theta_left_idx,
                                          NA),
                                   ifelse(self$reduced_model$theta_right_idx <=
                                            length(model$name_threshold),
                                          self$reduced_model$theta_right_idx,
                                          NA))],
                           ifelse(
                             self$reduced_model$theta_matrix_idx == 1,
                             self$reduced_model$theta_left_idx,
                             ifelse(
                               self$reduced_model$theta_matrix_idx == 2,
                               self$reduced_model$n_eta *
                                 (self$reduced_model$theta_left_idx - 1L) +
                                 self$reduced_model$theta_right_idx,
                               ifelse(self$reduced_model$theta_matrix_idx == 3,
                                      as.integer(
                                        self$reduced_model$n_eta *
                                          (self$reduced_model$theta_left_idx - 1L) +
                                          self$reduced_model$theta_right_idx
                                      ),
                                      as.integer(self$reduced_model$theta_right_idx))
                             )
                           ))
                  
                  self$reduced_model$n_theta_is_free <-
                    sum(self$reduced_model$theta_is_free)
                  self$reduced_model$n_theta_is_pen <-
                    sum(self$reduced_model$theta_is_pen)
                  
                  idx_row <- matrix(1:self$reduced_model$n_eta, 
                                self$reduced_model$n_eta, 
                                self$reduced_model$n_eta)
                  idx_col <- matrix(1:self$reduced_model$n_eta, 
                                    self$reduced_model$n_eta, 
                                    self$reduced_model$n_eta, 
                                    byrow = TRUE)
                  idx_matrix <- matrix(1:(self$reduced_model$n_eta^2), 
                                self$reduced_model$n_eta, 
                                self$reduced_model$n_eta, byrow = TRUE)
                  self$reduced_model$idx_vech <- 
                    matrix(1:(self$reduced_model$n_response^2), 
                           self$reduced_model$n_response,
                           self$reduced_model$n_response)[
                             lower.tri(matrix(0, 
                                              self$reduced_model$n_response, 
                                              self$reduced_model$n_response), 
                                       diag = TRUE)]
                  self$reduced_model$idx_nd_vech <- which(idx_row > idx_col)
                  self$reduced_model$idx_nd_tvech <- idx_matrix[idx_row > idx_col]
                })

## \code{$initialize_reduced_data()} initializes a reduced data. ##
lslxFitting$set("private",
                "initialize_reduced_data",
                function(data) {
                  self$reduced_data <-
                    list(
                      n_observation = integer(),
                      n_complete_observation = integer(),
                      n_missing_pattern = integer(),
                      n_response = length(data$numeric_variable) + 
                        length(data$ordered_variable),
                      n_numeric = length(data$numeric_variable),
                      n_ordered = length(data$ordered_variable),
                      n_group = length(data$level_group),
                      n_threshold = sum(data$nlevel_ordered - 1),
                      name_response = data$name_response,
                      name_numeric = data$name_numeric,
                      name_ordered = data$name_ordered,
                      sample_proportion = list(),
                      saturated_threshold = list(),
                      saturated_mean = list(),
                      saturated_cov = list(),
                      saturated_moment = list(),
                      saturated_moment_acov = list()
                    )
                  if (self$control$response) {
                    idc_subset <-
                      lapply(
                        data$response,
                        FUN = function(response_i) {
                          idc_subset_i <-
                            (row.names(response_i) %in% self$control$subset)
                          return(idc_subset_i)
                        }
                      )
                    idc_complete <-
                      lapply(
                        X = data$pattern,
                        FUN = function(pattern_i) {
                          idc_complete_i <-
                            apply(X = pattern_i,
                                  MARGIN = 1,
                                  FUN = prod)
                          return(as.logical(idc_complete_i))
                        }
                      )
                    if (self$control$missing_method %in% c("two_stage")) {
                      idc_use <-
                        mapply(
                          FUN = function(pattern_i,
                                         idc_subset_i) {
                            idc_use_i <-
                              ((rowSums(pattern_i) > 0) &
                                 idc_subset_i)
                            return(idc_use_i)
                          },
                          data$pattern,
                          idc_subset,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else if (self$control$missing_method == "listwise_deletion") {
                      idc_use <-
                        mapply(
                          FUN = function(idc_complete_i,
                                         idc_subset_i) {
                            idc_use_i <- idc_complete_i & idc_subset_i
                            return(idc_use_i)
                          },
                          idc_complete,
                          idc_subset,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else {
                    }
                    if (self$control$auxiliary &
                        self$control$missing_method == "two_stage") {
                      response <-
                        mapply(
                          FUN = function(response_i,
                                         auxiliary_i,
                                         idc_use_i) {
                            return(cbind(response_i[idc_use_i, , drop = FALSE],
                                         auxiliary_i[idc_use_i, , drop = FALSE]))
                          },
                          data$response,
                          data$auxiliary,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      pattern <-
                        mapply(
                          FUN = function(pattern_i,
                                         auxiliary_i,
                                         idc_use_i) {
                            return(cbind(pattern_i[idc_use_i, , drop = FALSE],!is.na(auxiliary_i[idc_use_i, , drop = FALSE])))
                          },
                          data$pattern,
                          data$auxiliary,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else {
                      response <-
                        mapply(
                          FUN = function(response_i,
                                         idc_use_i) {
                            return(response_i[idc_use_i, , drop = FALSE])
                          },
                          data$response,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      pattern <-
                        mapply(
                          FUN = function(pattern_i,
                                         idc_use_i) {
                            return(pattern_i[idc_use_i, , drop = FALSE])
                          },
                          data$pattern,
                          idc_use,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    }
                    weight <-
                      mapply(
                        FUN = function(weight_i,
                                       idc_use_i) {
                          weight_i <- weight_i[idc_use_i, ]
                          weight_i <- weight_i / sum(weight_i)
                          return(weight_i)
                        },
                        data$weight,
                        idc_use,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE
                      )
                    self$reduced_data$n_observation <-
                      sum(sapply(X = response, FUN = nrow))
                    self$reduced_data$n_complete_observation <-
                      sum(unlist(idc_complete) & unlist(idc_use))
                    self$reduced_data$sample_proportion <-
                      lapply(
                        X = lapply(X = response, FUN = nrow),
                        FUN = function(sample_size_i) {
                          sample_proportion_i <-
                            sample_size_i / self$reduced_data$n_observation
                          return(sample_proportion_i)
                        }
                      )
                    m_factor <-
                      lapply(
                        X = pattern,
                        FUN = function(pattern_i) {
                          as.factor(apply(
                            X = pattern_i,
                            MARGIN = 1,
                            FUN = function(pattern_ij) {
                              do.call(what = paste,
                                      args = as.list(c(pattern_ij, sep = "/")))
                            }
                          ))
                        }
                      )
                    self$reduced_data$n_missing_pattern <-
                      nlevels(unlist(m_factor))
                    
                    if (self$control$continuous) {
                      self$reduced_data$saturated_mean <-
                        lapply(
                          X = response,
                          FUN = function(response_i) {
                            as.matrix(colMeans(response_i, na.rm = TRUE))
                          }
                        )
                      self$reduced_data$saturated_cov <-
                        lapply(X = response,
                               FUN = cov,
                               use = "pairwise.complete.obs")
                      m_idx <-
                        lapply(
                          X = m_factor,
                          FUN = function(m_factor_i) {
                            m_idx_i <-
                              lapply(
                                X = strsplit(levels(m_factor_i),
                                             split = "/"),
                                FUN = function(x) {
                                  return(as.integer(which(as.logical(x))) - 1)
                                }
                              )
                            names(m_idx_i) <- levels(m_factor_i)
                            return(m_idx_i)
                          }
                        )
                      m2_idx <-
                        lapply(
                          X = m_factor,
                          FUN = function(m_factor_i) {
                            m2_idx_i <-
                              lapply(
                                X = strsplit(levels(m_factor_i),
                                             split = "/"),
                                FUN = function(x) {
                                  m2_idx_i <- tcrossprod(as.matrix(as.logical(x)))
                                  m2_idx_i <-
                                    as.logical(m2_idx_i[lower.tri(m2_idx_i, diag = TRUE)])
                                  return(as.integer(which(m2_idx_i)) - 1)
                                }
                              )
                            names(m2_idx_i) <- levels(m_factor_i)
                            return(m2_idx_i)
                          }
                        )
                      y_obs <-
                        mapply(
                          FUN = function(response_i,
                                         m_factor_i,
                                         m_idx_i) {
                            y_i <- split(response_i, m_factor_i)
                            mapply(
                              FUN = function(y_ij,
                                             m_idx_ij) {
                                y_obs_ij <- as.matrix(y_ij[, (m_idx_ij + 1) , drop = FALSE])
                                storage.mode(y_obs_ij) <- "double"
                                return(y_obs_ij)
                              },
                              y_i,
                              m_idx_i,
                              SIMPLIFY = FALSE,
                              USE.NAMES = TRUE
                            )
                          },
                          response,
                          m_factor,
                          m_idx,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      w <-
                        mapply(
                          FUN = function(weight_i,
                                         m_factor_i) {
                            split(weight_i, m_factor_i)
                          },
                          weight,
                          m_factor,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      compute_saturated_moment_cpp(
                        y_obs = y_obs,
                        w = w,
                        m_idx = m_idx,
                        saturated_mean = self$reduced_data$saturated_mean,
                        saturated_cov = self$reduced_data$saturated_cov,
                        iter_other_max = self$control$iter_other_max,
                        tol_other = self$control$tol_other
                      )
                      self$reduced_data$saturated_moment_acov <-
                        lapply(
                          X = self$reduced_data$saturated_cov,
                          FUN = function(saturated_cov_i) {
                            
                          }
                        )
                      compute_saturated_moment_acov_response_cpp(
                        y_obs = y_obs,
                        w = w,
                        m_idx = m_idx,
                        m2_idx = m2_idx,
                        saturated_mean = self$reduced_data$saturated_mean,
                        saturated_cov = self$reduced_data$saturated_cov,
                        saturated_moment_acov = self$reduced_data$saturated_moment_acov
                      )
                    } else {
                      if (self$control$missing_method == "two_stage") {
                        stop(
                          "When some response variable is ordered, 'two_stage' method for missingness is not available."
                        )
                      } else {
                        level_group <- as.list(names(response))
                        names(level_group) <- level_group
                        lav_data <-
                          mapply(
                            FUN = function(response_i,
                                           level_group_i) {
                              response_i$.group <- level_group_i
                              return(as.data.frame(response_i))
                            },
                            response,
                            level_group,
                            SIMPLIFY = FALSE
                          )
                        if (length(response) > 1) {
                          lav_data <- do.call(rbind.data.frame, lav_data)
                        } else {
                          lav_data <- lav_data[[1]]
                        }
                        lav_fit <- lavaan::lavCor(
                          lav_data,
                          meanstructure = TRUE,
                          group = ".group",
                          missing = "listwise",
                          se = "robust.huber.white",
                          output = "fit"
                        )
                        saturated_moment <- lavaan::fitted(lav_fit)
                        if (length(response) > 1) {
                          saturated_moment <- saturated_moment[names(response)]
                        } else {
                          saturated_moment <- list(saturated_moment)
                          names(saturated_moment) <- names(response)
                        }
                        self$reduced_data$saturated_threshold <-
                          lapply(
                            FUN = function(saturated_moment_i) {
                              response_for_threshold_i <- 
                                sapply(strsplit(x = names(saturated_moment_i$th[]), 
                                                split = "\\|"),
                                       FUN = function(x) {"["(x, 1)})
                              mapply(FUN = function(name_response_i) {
                                idc_response_for_threshold_i <- 
                                  response_for_threshold_i %in% name_response_i
                                return(saturated_moment_i$th[][idc_response_for_threshold_i])
                              },
                              self$reduced_data$name_response,
                              SIMPLIFY = FALSE)
                            },
                            saturated_moment
                          )
                        self$reduced_data$saturated_mean <-
                          lapply(
                            X = saturated_moment,
                            FUN = function(saturated_moment_i) {
                              return(as.matrix(saturated_moment_i$mean[]))
                            }
                          )
                        self$reduced_data$saturated_cov <-
                          lapply(
                            X = saturated_moment,
                            FUN = function(saturated_moment_i) {
                              return(as.matrix(saturated_moment_i$cov[]))
                            }
                          )
                        self$reduced_data$saturated_moment_acov <-
                          lav_fit@SampleStats@NACOV
                        if (length(response) > 1) {
                          names(self$reduced_data$saturated_moment_acov) <-
                            lav_fit@Data@group.label
                          self$reduced_data$saturated_moment_acov <-
                            self$reduced_data$saturated_moment_acov[names(self$reduced_data$saturated_cov)]
                        } else {
                          names(self$reduced_data$saturated_moment_acov) <-
                            names(self$reduced_data$saturated_cov)
                        }
                        self$reduced_data$saturated_moment_acov <-
                          mapply(
                            FUN = function(saturated_moment_acov_i,
                                           sample_proportion_i,
                                           n_observation_i) {
                              saturated_moment_acov_i <-
                                saturated_moment_acov_i /
                                (sample_proportion_i * n_observation_i)
                            },
                            self$reduced_data$saturated_moment_acov,
                            self$reduced_data$sample_proportion,
                            self$reduced_data$n_observation,
                            SIMPLIFY = FALSE
                          )
                      }
                    }
                  } else {
                    self$reduced_data$n_observation <-
                      sum(unlist(data$sample_size))
                    self$reduced_data$n_complete_observation <-
                      self$reduced_data$n_observation
                    self$reduced_data$sample_proportion <-
                      lapply(
                        X = data$sample_size,
                        FUN = function(sample_size_i) {
                          sample_proportion_i <-
                            sample_size_i / self$reduced_data$n_observation
                          return(sample_proportion_i)
                        }
                      )
                    self$reduced_data$n_missing_pattern <- 1
                    self$reduced_data$saturated_mean <-
                      lapply(
                        X = data$sample_mean,
                        FUN = function(sample_mean_i) {
                          sample_mean_i <- as.matrix(sample_mean_i)
                          return(sample_mean_i)
                        }
                      )
                    self$reduced_data$saturated_cov <-
                      lapply(
                        X = data$sample_cov,
                        FUN = function(sample_cov_i) {
                          diag(sample_cov_i) <-
                            diag(sample_cov_i) + self$control$ridge_cov
                          return(sample_cov_i)
                        }
                      )
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_cov,
                        FUN = function(saturated_cov_i) {
                          
                        }
                      )
                    compute_saturated_moment_acov_moment_cpp(
                      n_observation = self$reduced_data$n_observation,
                      sample_proportion = self$reduced_data$sample_proportion,
                      saturated_cov = self$reduced_data$saturated_cov,
                      saturated_moment_acov = self$reduced_data$saturated_moment_acov
                    )
                  }
                  if (self$control$auxiliary &
                      self$control$missing_method == "two_stage") {
                    y_name <-
                      c(colnames(data$response[[1]]), colnames(data$auxiliary[[1]]))
                  } else {
                    if (self$control$response) {
                      y_name <- colnames(data$response[[1]])
                    } else {
                      y_name <- colnames(data$sample_cov[[1]])
                    }
                  }
                  mean_name <- paste0(y_name, "<-1")
                  cov_name <- outer(
                    y_name,
                    y_name,
                    FUN = function(y_name_i, y_name_j) {
                      return(paste(y_name_i, y_name_j, sep = "<->"))
                    }
                  )
                  self$reduced_data$saturated_mean <-
                    lapply(
                      X = self$reduced_data$saturated_mean,
                      FUN = function(saturated_mean_i) {
                        rownames(saturated_mean_i) <- y_name
                        return(saturated_mean_i)
                      }
                    )
                  self$reduced_data$saturated_cov <-
                    lapply(
                      X = self$reduced_data$saturated_cov,
                      FUN = function(saturated_cov_i) {
                        diag(saturated_cov_i) <-
                          diag(saturated_cov_i) + self$control$ridge_cov
                        rownames(saturated_cov_i) <- y_name
                        colnames(saturated_cov_i) <- y_name
                        return(saturated_cov_i)
                      }
                    )
                  if (self$control$continuous) {
                    cov_name <-
                      cov_name[lower.tri(cov_name, diag = TRUE)]
                    self$reduced_data$saturated_moment <-
                      mapply(
                        FUN = function(saturated_mean_i,
                                       saturated_cov_i) {
                          saturated_mean_i <- c(saturated_mean_i)
                          names(saturated_mean_i) <- mean_name
                          saturated_cov_i <-
                            saturated_cov_i[lower.tri(saturated_cov_i, diag = TRUE)]
                          names(saturated_cov_i) <- cov_name
                          saturated_moment_i <- 
                            c(saturated_mean_i, saturated_cov_i)
                          return(saturated_moment_i)
                        },
                        self$reduced_data$saturated_mean,
                        self$reduced_data$saturated_cov,
                        SIMPLIFY = FALSE
                      )
                    moment_name <- c(mean_name, cov_name)
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_moment_acov,
                        FUN = function(saturated_moment_acov_i) {
                          rownames(saturated_moment_acov_i) <-
                            moment_name
                          colnames(saturated_moment_acov_i) <-
                            moment_name
                          return(saturated_moment_acov_i)
                        }
                      )
                    if (self$control$auxiliary &
                        self$control$missing_method == "two_stage") {
                      y_name <- colnames(data$response[[1]])
                      mean_name <- paste0(y_name, "<-1")
                      cov_name <- outer(
                        y_name,
                        y_name,
                        FUN = function(y_name_i, y_name_j) {
                          return(paste(y_name_i, y_name_j, sep = "<->"))
                        }
                      )
                      cov_name <-
                        cov_name[lower.tri(cov_name, diag = TRUE)]
                      moment_name <- c(mean_name, cov_name)
                      self$reduced_data$saturated_mean <-
                        lapply(
                          X = self$reduced_data$saturated_mean,
                          FUN = function(saturated_mean_i) {
                            return(saturated_mean_i[y_name, , drop = FALSE])
                          }
                        )
                      self$reduced_data$saturated_cov <-
                        lapply(
                          X = self$reduced_data$saturated_cov,
                          FUN = function(saturated_cov_i) {
                            return(saturated_cov_i[y_name, y_name, drop = FALSE])
                          }
                        )
                      self$reduced_data$saturated_moment_acov <-
                        lapply(
                          X = self$reduced_data$saturated_moment_acov,
                          FUN = function(saturated_moment_acov_i) {
                            return(saturated_moment_acov_i[moment_name,
                                                           moment_name,
                                                           drop = FALSE])
                          }
                        )
                    } 
                  } else {
                    cov_name <-
                      cov_name[lower.tri(cov_name, diag = FALSE)]
                    self$reduced_data$saturated_moment <-
                      mapply(
                        FUN = function(saturated_threshold_i,
                                       saturated_mean_i,
                                       saturated_cov_i) {
                          rownames(saturated_mean_i) <- mean_name
                          if (length(data$numeric_variable) > 0) {
                            saturated_threshold_i <- 
                              split(
                                c(unlist(unname(saturated_threshold_i)), 
                                  - saturated_mean_i[match(data$numeric_variable,
                                                         data$name_response), 1]),
                                unlist(lapply(
                                  X = strsplit(names(c(unlist(unname(saturated_threshold_i)), 
                                                       saturated_mean_i[match(data$numeric_variable,
                                                                              data$name_response), 1])), "\\||<-"),
                                  FUN = function(x) {
                                    "["(x, 1)
                                  }
                                ))
                              )[data$name_response]
                            saturated_threshold_i <- unlist(unname(saturated_threshold_i))
                            saturated_var_i <-
                              diag(saturated_cov_i)[match(data$numeric_variable,
                                                          data$name_response)]
                            names(saturated_var_i) <-
                              paste0(data$numeric_variable,
                                     "<->",
                                     data$numeric_variable)
                          } else {
                            saturated_threshold_i <- unlist(unname(saturated_threshold_i))
                            saturated_var_i <- numeric(0)
                          }
                          saturated_cov_i <-
                            saturated_cov_i[lower.tri(saturated_cov_i, diag = FALSE)]
                          names(saturated_cov_i) <- cov_name
                          saturated_moment_i <- 
                            c(saturated_threshold_i, saturated_var_i, saturated_cov_i)
                          return(saturated_moment_i)
                        },
                        self$reduced_data$saturated_threshold,
                        self$reduced_data$saturated_mean,
                        self$reduced_data$saturated_cov,
                        SIMPLIFY = FALSE
                      )
                    moment_name <- 
                      names(self$reduced_data$saturated_moment[[1]])
                    self$reduced_data$saturated_moment_acov <-
                      lapply(
                        X = self$reduced_data$saturated_moment_acov,
                        FUN = function(saturated_moment_acov_i) {
                          rownames(saturated_moment_acov_i) <- 
                            moment_name
                          colnames(saturated_moment_acov_i) <- 
                            moment_name
                          return(saturated_moment_acov_i)
                        }
                      )
                  }
                  
                })


## \code{$initialize_weight_matrix()} initializes weight matrix. ##
lslxFitting$set("private",
                "initialize_weight_matrix",
                function() {
                  if (self$control$loss %in% c("uls", "dwls", "wls")) {
                    if (is.null(self$control$weight_matrix)) {
                      self$control$weight_matrix <-
                        mapply(
                          FUN = function(saturated_moment_acov_i,
                                         sample_proportion_i) {
                            if (self$control$loss == "uls") {
                              weight_matrix_i <-
                                sample_proportion_i * diag(dim(saturated_moment_acov_i)[1])
                            } else if (self$control$loss == "dwls") {
                              weight_matrix_i <-
                                diag(1 / diag(saturated_moment_acov_i)) / self$reduced_data$n_observation
                            } else if (self$control$loss == "wls") {
                              weight_matrix_i <-
                                solve(saturated_moment_acov_i) / self$reduced_data$n_observation
                            } else {
                            }
                            return(weight_matrix_i)
                          },
                          self$reduced_data$saturated_moment_acov,
                          self$reduced_data$sample_proportion,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                    } else {
                      if (!all(dim(self$control$weight_matrix[[1]]) == 
                               c(self$reduced_model$n_moment, 
                                 self$reduced_model$n_moment))) {
                        stop(
                          "The dimension of some element in 'weight_matrix' is not correct.",
                          "\n  The correct dimension is ",
                          n_moment,
                          " by ",
                          n_moment,
                          "."
                        )
                      }
                    }
                    self$control$weight_matrix <-
                      lapply(
                        X = self$control$weight_matrix,
                        FUN = function(weight_matrix_i) {
                          rownames(weight_matrix_i) <-
                            rownames(self$reduced_data$saturated_moment_acov[[1]])
                          colnames(weight_matrix_i) <-
                            colnames(self$reduced_data$saturated_moment_acov[[1]])
                          return(weight_matrix_i)
                        }
                      )
                  }
                })

## \code{$initialize_supplied_result()} initializes a supplied result. ##
lslxFitting$set("private",
                "initialize_supplied_result",
                function() {
                  self$supplied_result <- list()
                  private$compute_baseline_model()
                  private$compute_saturated_model()
                  private$compute_fitted_start()
                })

## \code{$compute_fitted_start()} computes starting value. ##
lslxFitting$set("private",
                "compute_fitted_start",
                function() {
                  self$supplied_result$fitted_start <- 
                    self$reduced_model$theta_start
                  if (!self$control$continuous) {
                    saturated_threshold_pool <-
                      mapply(
                        FUN = function(sample_proportion_i,
                                       saturated_threshold_i) {
                          mapply(
                            FUN = function(saturated_threshold_ij) {
                              return(sample_proportion_i * saturated_threshold_ij)  
                            },
                            saturated_threshold_i,
                            SIMPLIFY = FALSE,
                            USE.NAMES = TRUE
                          )
                        },
                        self$reduced_data$sample_proportion,
                        self$reduced_data$saturated_threshold,
                        SIMPLIFY = FALSE,
                        USE.NAMES = TRUE
                      )
                    if (self$reduced_data$n_group > 1) {
                      for (i in 2:self$reduced_data$n_group) {
                        saturated_threshold_pool[[1]] <-
                          Map("+", saturated_threshold_pool[[1]], 
                              saturated_threshold_pool[[i]])
                      }
                    } 
                    saturated_threshold_pool <- saturated_threshold_pool[[1]]
                    if (0 %in% unique(self$reduced_model$theta_group_idx)) {
                      idc_0 <-
                        (self$reduced_model$theta_matrix_idx == 0) &
                        (self$reduced_model$theta_group_idx == 0)
                      
                      self$supplied_result$fitted_start[idc_0] <-
                        ifelse(!is.na(self$reduced_model$theta_start[idc_0]),
                               self$reduced_model$theta_start[idc_0],
                               apply(cbind(
                                 self$reduced_model$theta_left_idx[idc_0],
                                 self$reduced_model$theta_right_idx[idc_0]
                               ),
                               MARGIN = 1,
                               function(idx) {
                                 saturated_threshold_pool[[idx[1]]][idx[2]]
                               }))
                      for (i in setdiff(unique(self$reduced_model$theta_group_idx), 0)) {
                        idc_i <-
                          (self$reduced_model$theta_matrix_idx == 0) &
                          (self$reduced_model$theta_group_idx == i)
                        self$supplied_result$fitted_start[idc_i] <-
                          ifelse(
                            is.na(self$reduced_model$theta_start[idc_i]),
                            0,
                            self$reduced_model$theta_start[idc_i]
                          )
                      }
                    } else {
                      for (i in seq_len(self$reduced_model$n_group)) {
                        idc_i <-
                          (self$reduced_model$theta_matrix_idx == 0) &
                          (self$reduced_model$theta_group_idx == i)
                        self$supplied_result$fitted_start[idc_i] <-
                          ifelse(!is.na(self$reduced_model$theta_start[idc_i]),
                                 self$reduced_model$theta_start[idc_i],
                                 apply(cbind(
                                   self$reduced_model$theta_left_idx[idc_i],
                                   self$reduced_model$theta_right_idx[idc_i]
                                 ),
                                 MARGIN = 1,
                                 function(idx) {
                                   saturated_threshold_pool[[idx[1]]][idx[2]]
                                 }))
                      }
                    }
                  }
                  
                  if (self$control$start_method == "mh") {
                    saturated_cov_pool <-
                      Reduce(
                        "+",
                        mapply(
                          FUN = function(sample_proportion_i,
                                         saturated_cov_i) {
                            return(sample_proportion_i * saturated_cov_i)
                          },
                          self$reduced_data$sample_proportion,
                          self$reduced_data$saturated_cov,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      )
                    saturated_mean_pool <-
                      Reduce(
                        "+",
                        mapply(
                          FUN = function(sample_proportion_i,
                                         saturated_mean_i) {
                            return(sample_proportion_i * saturated_mean_i)
                          },
                          self$reduced_data$sample_proportion,
                          self$reduced_data$saturated_mean,
                          SIMPLIFY = FALSE,
                          USE.NAMES = TRUE
                        )
                      )
                    if (self$control$model_class == "efa") {
                      idc_beta <-
                        (self$reduced_model$theta_matrix_idx == 2 &
                           self$reduced_model$theta_is_pen)
                    } else {
                      idc_beta <-
                        (self$reduced_model$theta_matrix_idx == 2 &
                           self$reduced_model$theta_is_free) |
                        (
                          self$reduced_model$theta_matrix_idx == 2 &
                            (
                              !self$reduced_model$theta_is_free &
                                !self$reduced_model$theta_is_pen
                            ) &
                            (self$reduced_model$theta_start != 0)
                        ) |
                        (
                          self$reduced_model$theta_matrix_idx == 2 &
                            self$reduced_model$theta_is_pen &
                            !(
                              self$reduced_model$theta_left_idx <=
                                self$reduced_model$n_response &
                                self$reduced_model$theta_right_idx >
                                self$reduced_model$n_response
                            )
                        )
                    }
                    idx_beta <-
                      strsplit(x = unique(
                        paste0(
                          self$reduced_model$theta_left_idx[idc_beta],
                          ",",
                          self$reduced_model$theta_right_idx[idc_beta]
                        )
                      ),
                      split = ",")
                    theta_left_idx_beta <-
                      sapply(
                        X = idx_beta,
                        FUN = function(idx_beta_i) {
                          as.integer(idx_beta_i[1])
                        }
                      )
                    theta_right_idx_beta <-
                      sapply(
                        X = idx_beta,
                        FUN = function(idx_beta_i) {
                          as.integer(idx_beta_i[2])
                        }
                      )
                    if (self$reduced_model$n_factor > 0) {
                      cov_eta <-
                        matrix(0,
                               self$reduced_model$n_eta,
                               self$reduced_model$n_eta)
                      cor_pool <- cov2cor(saturated_cov_pool)
                      cov_eta[1:self$reduced_model$n_response,
                              1:self$reduced_model$n_response] <-
                        cor_pool
                      if (self$control$model_class == "efa") {
                        cor_pool_eigen <- eigen(cor_pool)
                        cov_eta_yf <-
                          cor_pool_eigen$vectors[, 1:self$reduced_model$n_factor, drop = FALSE] %*%
                          diag(sqrt(cor_pool_eigen$values[1:self$reduced_model$n_factor]))
                        cov_eta_yf <-
                          promax(cov_eta_yf, m = 4)$loadings[]
                        cov_eta_yf[abs(cov_eta_yf) < 0.2] <- 0
                        cov_eta[1:self$reduced_model$n_response,
                                (self$reduced_model$n_response + 1):self$reduced_model$n_eta] <-
                          cov_eta_yf
                        cov_eta[(self$reduced_model$n_response + 1):self$reduced_model$n_eta,
                                1:self$reduced_model$n_response] <-
                          t(cov_eta_yf)
                        cov_eta[(self$reduced_model$n_response + 1):self$reduced_model$n_eta,
                                (self$reduced_model$n_response + 1):self$reduced_model$n_eta] <-
                          diag(1, self$reduced_model$n_factor)
                      } else {
                        for (i in (self$reduced_model$n_response + 1):(self$reduced_model$n_eta)) {
                          theta_left_idx_beta_i <-
                            theta_left_idx_beta[theta_right_idx_beta == i &
                                                  theta_left_idx_beta <= self$reduced_model$n_response]
                          cov_sum_i <-
                            sum(cov_eta[theta_left_idx_beta_i,
                                        theta_left_idx_beta_i])
                          for (j in seq_len(i)) {
                            if (j <= self$reduced_model$n_response) {
                              cov_eta_ij <-
                                sum(cov_eta[j, theta_left_idx_beta_i]) / sqrt(abs(cov_sum_i))
                              cov_eta[i, j] <- cov_eta_ij
                              cov_eta[j, i] <- cov_eta_ij
                            } else if (j < i) {
                              theta_left_idx_beta_j <-
                                theta_left_idx_beta[theta_right_idx_beta == j &
                                                      theta_left_idx_beta <= self$reduced_model$n_response]
                              cov_sum_j <-
                                sum(cov_eta[theta_left_idx_beta_j,
                                            theta_left_idx_beta_j])
                              cov_eta_ij <-
                                sum(cov_eta[theta_left_idx_beta_j, theta_left_idx_beta_i]) /
                                sqrt(abs(cov_sum_i) * abs(cov_sum_j))
                              cov_eta[i, j] <- cov_eta_ij
                              cov_eta[j, i] <- cov_eta_ij
                            } else {
                              cov_eta[j, i] <- 1
                            }
                          }
                        }
                      }
                      cov_eta <-
                        diag(c(sqrt(diag(
                          saturated_cov_pool
                        )),
                        rep(1, self$reduced_model$n_factor))) %*%
                        cov_eta %*%
                        diag(c(sqrt(diag(
                          saturated_cov_pool
                        )),
                        rep(1, self$reduced_model$n_factor)))
                    } else {
                      cov_eta <- saturated_cov_pool
                    }
                    beta_start <-
                      matrix(0,
                             self$reduced_model$n_eta,
                             self$reduced_model$n_eta)
                    
                    for (i in seq_len(self$reduced_model$n_eta)) {
                      theta_right_idx_beta_i <-
                        theta_right_idx_beta[theta_left_idx_beta == i]
                      if (length(theta_right_idx_beta_i) > 0) {
                        beta_start_i <-
                          solve(cov_eta[theta_right_idx_beta_i,
                                        theta_right_idx_beta_i,
                                        drop = FALSE]) %*%
                          cov_eta[theta_right_idx_beta_i, i, drop = FALSE]
                        beta_start[i, theta_right_idx_beta_i] <-
                          beta_start_i
                      }
                    }
                    idc_phi <-
                      (self$reduced_model$theta_matrix_idx == 3 &
                         self$reduced_model$theta_is_free) |
                      (
                        self$reduced_model$theta_matrix_idx == 3 &
                          (
                            !self$reduced_model$theta_is_free &
                              !self$reduced_model$theta_is_pen
                          ) &
                          (self$reduced_model$theta_start != 0)
                      )
                    idc_phi[is.na(idc_phi)] <- TRUE
                    idx_phi <-
                      strsplit(x = unique(
                        paste0(
                          self$reduced_model$theta_left_idx[idc_phi],
                          ",",
                          self$reduced_model$theta_right_idx[idc_phi]
                        )
                      ),
                      split = ",")
                    theta_left_idx_phi <-
                      sapply(
                        X = idx_phi,
                        FUN = function(idx_phi_i) {
                          as.integer(idx_phi_i[1])
                        }
                      )
                    theta_right_idx_phi <-
                      sapply(
                        X = idx_phi,
                        FUN = function(idx_phi_i) {
                          as.integer(idx_phi_i[2])
                        }
                      )
                    phi_start <-
                      matrix(0,
                             self$reduced_model$n_eta,
                             self$reduced_model$n_eta)
                    phi_saturated <-
                      (diag(1, self$reduced_model$n_eta) - beta_start) %*%
                      cov_eta %*% t((diag(1, self$reduced_model$n_eta) - beta_start))
                    phi_start[cbind(theta_left_idx_phi, theta_right_idx_phi)] <-
                      phi_saturated[cbind(theta_left_idx_phi, theta_right_idx_phi)]
                    phi_start[cbind(theta_right_idx_phi, theta_left_idx_phi)] <-
                      phi_saturated[cbind(theta_right_idx_phi, theta_left_idx_phi)]
                    idc_alpha <-
                      (self$reduced_model$theta_matrix_idx == 1 &
                         self$reduced_model$theta_is_free) |
                      (
                        self$reduced_model$theta_matrix_idx == 1 &
                          (
                            !self$reduced_model$theta_is_free &
                              !self$reduced_model$theta_is_pen
                          ) &
                          (self$reduced_model$theta_start != 0)
                      )
                    theta_left_idx_alpha <-
                      unique(self$reduced_model$theta_left_idx[idc_alpha])
                    if (self$reduced_model$n_factor > 0) {
                      mean_eta <-
                        matrix(0, self$reduced_model$n_eta, 1)
                      mean_eta[1:self$reduced_model$n_response, 1] <-
                        saturated_mean_pool
                      for (i in (self$reduced_model$n_response + 1):(self$reduced_model$n_eta)) {
                        if (i %in% theta_left_idx_alpha) {
                          theta_left_idx_beta_i <-
                            theta_left_idx_beta[theta_right_idx_beta == i]
                          mean_eta[i, 1] <-
                            solve(t(beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                                    beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                            t(beta_start[theta_left_idx_beta_i, i, drop = FALSE]) %*%
                            mean_eta[theta_left_idx_beta_i, 1, drop = FALSE]
                        } else {
                          mean_eta[i, 1] <- 0
                        }
                      }
                    } else {
                      mean_eta <- saturated_mean_pool
                    }
                    alpha_saturated <-
                      (diag(1, self$reduced_model$n_eta) - beta_start) %*% mean_eta
                    alpha_start <-
                      matrix(0, self$reduced_model$n_eta, 1)
                    alpha_start[theta_left_idx_alpha, 1] <-
                      alpha_saturated[theta_left_idx_alpha, 1]
                    
                    coefficient_matrix_start <-
                      list(
                        alpha_start = alpha_start,
                        beta_start = beta_start,
                        phi_start = phi_start
                      )
                    for (i in 1:3) {
                      if (0 %in% unique(self$reduced_model$theta_group_idx)) {
                        idc_i0 <-
                          (self$reduced_model$theta_matrix_idx == i) &
                          (self$reduced_model$theta_group_idx == 0)
                        
                        self$supplied_result$fitted_start[idc_i0] <-
                          ifelse(
                            is.na(self$reduced_model$theta_start[idc_i0]),
                            coefficient_matrix_start[[i]][cbind(self$reduced_model$theta_left_idx[idc_i0],
                                                                self$reduced_model$theta_right_idx[idc_i0])],
                            self$reduced_model$theta_start[idc_i0]
                          )
                        for (j in setdiff(unique(self$reduced_model$theta_group_idx), 0)) {
                          idc_ij <-
                            (self$reduced_model$theta_matrix_idx == i) &
                            (self$reduced_model$theta_group_idx == j)
                          self$supplied_result$fitted_start[idc_ij] <-
                            ifelse(
                              is.na(self$reduced_model$theta_start[idc_ij]),
                              0,
                              self$reduced_model$theta_start[idc_ij]
                            )
                        }
                      } else {
                        for (j in seq_len(self$reduced_model$n_group)) {
                          idc_ij <-
                            (self$reduced_model$theta_matrix_idx == i) &
                            (self$reduced_model$theta_group_idx == j)
                          self$supplied_result$fitted_start[idc_ij] <-
                            ifelse(
                              is.na(self$reduced_model$theta_start[idc_ij]),
                              coefficient_matrix_start[[i]][cbind(
                                self$reduced_model$theta_left_idx[idc_ij],
                                self$reduced_model$theta_right_idx[idc_ij]
                              )],
                              self$reduced_model$theta_start[idc_ij]
                            )
                        }
                      }
                    }
                    self$supplied_result$alpha_start <- alpha_start
                    self$supplied_result$beta_start <- beta_start
                    self$supplied_result$phi_start <- phi_start
                  } else if (self$control$start_method == "heuristic") {
                    self$supplied_result$fitted_start <-
                      ifelse(
                        !is.na(self$reduced_model$theta_start),
                        self$reduced_model$theta_start,
                        ifelse(
                          self$reduced_model$theta_matrix_idx == 2,
                          ifelse(self$reduced_model$theta_is_free,
                                 0.5, 0),
                          ifelse((self$reduced_model$theta_matrix_idx == 3) &
                                   (
                                     self$reduced_model$theta_left_idx ==
                                       self$reduced_model$theta_right_idx
                                   ),
                                 .1,
                                 0
                          )
                        )
                      )
                    if (any(self$reduced_model$theta_group_idx == 0)) {
                      self$supplied_result$fitted_start <-
                        ifelse(
                          self$reduced_model$theta_group_idx == 0,
                          self$supplied_result$fitted_start,
                          0
                        )
                    }
                  } else {
                    self$supplied_result$fitted_start <- self$reduced_model$theta_start
                  }
                  self$supplied_result$fitted_start[
                    is.na(self$supplied_result$fitted_start)] <- 0
                  names(self$supplied_result$fitted_start) <-
                    self$reduced_model$theta_name
                })

## \code{$compute_baseline_model()} computes baseline model. ##
lslxFitting$set("private",
                "compute_baseline_model",
                function() {
                  if (self$control$continuous) {
                      loss_value =
                        ifelse(self$control$loss == "ml",
                               sum(
                                 mapply(
                                   FUN = function(saturated_mean_i,
                                                  saturated_cov_i,
                                                  sample_proportion_i) {
                                     baseline_loss_value_i <-
                                       sample_proportion_i *
                                       (sum(diag(saturated_cov_i %*%
                                                   diag(
                                                     1 / diag(saturated_cov_i)
                                                   ))) -
                                          log(det(saturated_cov_i %*% diag(
                                            1 / diag(saturated_cov_i)
                                          ))) -
                                          self$reduced_model$n_response)
                                     return(baseline_loss_value_i)
                                   },
                                   self$reduced_data$saturated_mean,
                                   self$reduced_data$saturated_cov,
                                   self$reduced_data$sample_proportion,
                                   SIMPLIFY = TRUE
                                 )
                               ),
                               sum(
                                 mapply(
                                   FUN = function(saturated_mean_i,
                                                  saturated_cov_i,
                                                  weight_matrix_i) {
                                     model_residual_i <-
                                       as.matrix(c(
                                         saturated_mean_i - saturated_mean_i,
                                         diag(diag(saturated_cov_i))[lower.tri(saturated_cov_i, diag = TRUE)] -
                                           saturated_cov_i[lower.tri(saturated_cov_i, diag = TRUE)]
                                       ))
                                     baseline_loss_value_i <-
                                       t(model_residual_i) %*% weight_matrix_i %*% model_residual_i
                                     return(baseline_loss_value_i)
                                   },
                                   self$reduced_data$saturated_mean,
                                   self$reduced_data$saturated_cov,
                                   self$control$weight_matrix,
                                   SIMPLIFY = TRUE
                                 )
                               ))
                      n_nonzero_coefficient <-
                        self$reduced_model$n_group *
                        (2L * self$reduced_model$n_response)
                      degrees_of_freedom <-
                        self$reduced_model$n_group *
                        (
                          self$reduced_model$n_moment -
                            2L * self$reduced_model$n_response
                        )
                  } else {
                      loss_value <-
                        sum(
                          mapply(
                            FUN = function(saturated_moment_i,
                                           weight_matrix_i) {
                              model_residual_i <-
                                as.matrix(saturated_moment_i)
                              model_residual_i[
                                ((self$reduced_data$n_threshold + 
                                    length(self$reduced_data$n_numeric) * 2) + 
                                   1):nrow(model_residual_i), ] <- 0
                              baseline_loss_value_i <-
                                t(model_residual_i) %*% weight_matrix_i %*% model_residual_i
                              return(baseline_loss_value_i)
                            },
                            self$reduced_data$saturated_moment,
                            self$control$weight_matrix,
                            SIMPLIFY = TRUE
                          )
                        )
                      n_nonzero_coefficient <-
                        self$reduced_model$n_group *
                        (self$reduced_data$n_threshold + 
                           length(self$reduced_data$n_numeric) * 2L)
                      degrees_of_freedom <-
                        self$reduced_model$n_group *
                        (
                          self$reduced_model$n_moment -
                            (self$reduced_data$n_threshold + 
                               length(self$reduced_data$n_numeric) * 2L)
                        )
                  }
                  self$supplied_result$baseline_model <- 
                    c(loss_value = loss_value, 
                      n_nonzero_coefficient = n_nonzero_coefficient,
                      degrees_of_freedom = degrees_of_freedom)
                })

## \code{$compute_saturated_model()} computes saturated model. ##
lslxFitting$set("private",
                "compute_saturated_model",
                function() {
                  self$supplied_result$saturated_model <- c(
                    loss_value = 0,
                    n_nonzero_coefficient = self$reduced_model$n_group *
                      self$reduced_model$n_moment,
                    degrees_of_freedom = 0
                  )
                })


## \code{$initialize_grid()} initializes lambda_grid and delta_grid by default. ##
lslxFitting$set("private",
                "initialize_grid",
                function() {
                  if (any(self$reduced_model$eta_is_exogenous)) {
                    eta_is_exogenous <- self$reduced_model$eta_is_exogenous
                  } else {
                    eta_is_exogenous <- rep(TRUE, self$reduced_model$n_eta)
                  }
                  if (any(self$reduced_model$eta_is_endogenous)) {
                    eta_is_endogenous <- self$reduced_model$eta_is_endogenous
                  } else {
                    eta_is_endogenous <- rep(TRUE, self$reduced_model$n_eta)
                  }
                  if (self$control$regularizer) {
                    if (self$control$lambda_grid[[1]] == "default") {
                      if (self$control$start_method == "mh") {
                        beta_start_inv <-
                          solve(diag(self$reduced_model$n_eta) - self$supplied_result$beta_start)
                        sigma_eta_start <-
                          beta_start_inv %*% self$supplied_result$phi_start %*% t(beta_start_inv)
                        lambda_max <-
                          sqrt(quantile(diag(sigma_eta_start)[eta_is_exogenous], 0.6)) /
                          quantile(diag(self$supplied_result$phi_start)[eta_is_endogenous] /
                                     sqrt(diag(sigma_eta_start)[eta_is_endogenous]),
                                   0.4) *
                          self$control$threshold_value
                      } else if (self$control$start_method == "heuristic") {
                        saturated_var <- diag(do.call("+", self$reduced_data$saturated_cov))
                        lambda_max <-
                          1 / quantile((saturated_var * 0.3) / sqrt(saturated_var), 0.3) *
                          self$control$threshold_value
                      } else {}
                      lambda_min <-
                        (log(self$reduced_data$n_observation) / self$reduced_data$n_observation)
                      if (self$control$double_regularizer) {
                        self$control$lambda_grid <- 
                          list(exp(seq(
                            log(lambda_max),
                            log(lambda_min),
                            length.out = self$control$lambda_length
                          )),
                          exp(seq(
                            log(lambda_max),
                            log(lambda_min),
                            length.out = self$control$lambda_length
                          )))
                      } else {
                        self$control$lambda_grid <- 
                          list(exp(seq(
                            log(lambda_max),
                            log(lambda_min),
                            length.out = self$control$lambda_length
                          )),
                          0)
                      }
                    } else {
                      if (self$control$double_regularizer) {
                        self$control$lambda_grid <- 
                          list(self$control$lambda_grid,
                               self$control$lambda_grid)
                      } else {
                        self$control$lambda_grid <- 
                          list(self$control$lambda_grid, 0)
                      }
                    }
                    if (self$control$delta_grid[[1]] == "default") {
                      if (self$control$penalty_method == "lasso") {
                        self$control$delta_grid <- list(1, 1)
                      } else if (self$control$penalty_method == "ridge") {
                        self$control$delta_grid <- list(0, 0)
                      } else if (self$control$penalty_method == "elastic_net") {
                        self$control$delta_grid <-
                          seq(0, 1, length.out = (self$control$delta_length + 2))[2:(self$control$delta_length + 1)]
                        if (self$control$double_regularizer) {
                          self$control$delta_grid <-
                            list(self$control$delta_grid,
                                 self$control$delta_grid) 
                        } else {
                          self$control$delta_grid <-
                            list(self$control$delta_grid,
                                 0.5) 
                        }
                      } else if (self$control$penalty_method %in% c("mcp")) {
                        if (self$control$start_method == "mh") {
                          beta_start_inv <-
                            solve(diag(self$reduced_model$n_eta) - self$supplied_result$beta_start)
                          sigma_eta_start <-
                            beta_start_inv %*% self$supplied_result$phi_start %*% t(beta_start_inv)
                          delta_min <-
                            2 * max(diag(sigma_eta_start)[eta_is_endogenous]) /
                            min(diag(sigma_eta_start)[eta_is_exogenous])
                        } else if (self$control$start_method == "heuristic") {
                          saturated_var <- diag(do.call("+", self$reduced_data$saturated_cov))
                          delta_min <- 2 * max(saturated_var, 1) / 1
                        } else {
                          
                        }
                        if (self$control$delta_length == 1) {
                          self$control$delta_grid <- 
                            list(Inf, Inf)
                        } else {
                          if (self$control$double_regularizer) {
                            self$control$delta_grid <-
                              list(c(delta_min * (1:(self$control$delta_length - 1)), Inf),
                                   c(delta_min * (1:(self$control$delta_length - 1)), Inf))
                          } else {
                            self$control$delta_grid <-
                              list(c(delta_min * (1:(self$control$delta_length - 1)), Inf),
                                   Inf) 
                          }
                        }
                      } else {}
                    } else {
                      if (self$control$double_regularizer) {
                        self$control$delta_grid <- 
                          list(self$control$delta_grid,
                               self$control$delta_grid)
                      } else {
                        self$control$delta_grid <- 
                          list(self$control$delta_grid,
                               ifelse(self$control$penalty_method == "mcp", Inf, 0))
                      }
                    }
                    if (self$control$lambda_direction == "default") {
                      if (min(self$control$lambda_grid[[1]]) == 0) {
                        self$control$lambda_direction <- "decrease"
                      } else {
                        self$control$lambda_direction <- "increase"
                      }
                    }
                    if (self$control$lambda_direction == "decrease") {
                      self$control$lambda_grid <-
                        lapply(X = self$control$lambda_grid,
                               FUN = function(x) {
                                 sort(x, decreasing = TRUE)
                               })
                    } else if (self$control$lambda_direction == "increase") {
                      self$control$lambda_grid <-
                        lapply(X = self$control$lambda_grid,
                               FUN = function(x) {
                                 sort(x, decreasing = FALSE)
                               })
                    } else {}
                    if (self$control$penalty_method == "elastic_net") {
                      self$control$delta_grid <-
                        lapply(X = self$control$delta_grid,
                               FUN = function(x) {
                                 sort(x, decreasing = FALSE)
                               })
                    } else {
                      self$control$delta_grid <-
                        lapply(X = self$control$delta_grid,
                               FUN = function(x) {
                                 sort(x, decreasing = TRUE)
                               })
                    }
                  } else {
                    self$control$lambda_grid <- list(0, 0)
                    self$control$delta_grid <- list(ifelse(self$control$penalty_method == "mcp", Inf, 0), 
                                                    ifelse(self$control$penalty_method == "mcp", Inf, 0))
                  }
                  names(self$control$lambda_grid) <- c("lambda_1st_grid", "lambda_2nd_grid")
                  names(self$control$delta_grid) <- c("delta_1st_grid", "delta_2nd_grid")
                  if (self$control$searcher) {
                    if (self$control$step_grid[[1]] == "default") {
                      self$control$step_grid <-
                        0:self$reduced_model$n_theta_is_pen
                    } else {
                      self$control$step_grid <-
                        0:min(self$reduced_model$n_theta_is_pen,
                              max(round(self$control$step_grid)))
                    }
                  } else {
                    self$control$step_grid = 0
                  }
                })


## \code{$initialize_fitted_result()} initializes a fitted result. ##
lslxFitting$set("private",
                "initialize_fitted_result",
                function() {
                  self$fitted_result <- list()
                  if (self$control$regularizer) {
                    length_fitted_result <-
                      (length(self$control$lambda_grid[[1]]) *
                         length(self$control$delta_grid[[1]])) *
                      (length(self$control$lambda_grid[[2]]) *
                         length(self$control$delta_grid[[2]]))
                  } else {
                    length_fitted_result <-
                      length(self$control$step_grid)
                  }
                  self$fitted_result$numerical_condition <-
                    vector(mode = "list",
                           length = length_fitted_result)
                  self$fitted_result$information_criterion <-
                    vector(mode = "list",
                           length = length_fitted_result)
                  self$fitted_result$fit_index <-
                    vector(mode = "list",
                           length = length_fitted_result)
                  self$fitted_result$cv_error <-
                    vector(mode = "list",
                           length = length_fitted_result)
                  self$fitted_result$coefficient <-
                    vector(mode = "list",
                           length = length_fitted_result)
                })
