## \code{$new()} initialzes a new object of \code{prelslx} R6 class.##
prelslx$set("public",
            "initialize",
            function(model,
                     data,
                     numeric_variable,
                     ordered_variable,
                     weight_variable,
                     auxiliary_variable,
                     group_variable,
                     reference_group,
                     sample_cov,
                     sample_mean,
                     sample_size,
                     sample_moment_acov,
                     verbose = TRUE) {
              if (missing(data) & missing(sample_cov)) {
                stop("Argument 'data' and 'sample_cov' cannot be both empty.")
              } else {
                if (!missing(data)) {
                  if (is.matrix(data) & is.numeric(data)) {
                    data <- as.data.frame(data)
                  }
                  if (!is.data.frame(data)) {
                    stop("Argument 'data' is not a 'data.frame'.")
                  }
                  if (!missing(numeric_variable) &
                      !missing(ordered_variable)) {
                    if (length(intersect(numeric_variable, ordered_variable)) == 0) {
                      stop("Arguments 'numeric_variable' and 'ordered_variable' share the same variables.")
                    }
                  }
                  if (!missing(numeric_variable)) {
                    if (!is.character(numeric_variable)) {
                      stop("Argument 'numeric_variable' is not a 'character'.")
                    }
                    if (!all(numeric_variable %in% colnames(data))) {
                      stop("Some element in argument 'numeric_variable' is not recognized.")
                    } else {
                      data[, numeric_variable] <- 
                        lapply(X = numeric_variable, 
                               FUN = function(numeric_variable_i) {
                                 return(as.numeric(getElement(data, numeric_variable_i)))
                               })
                    }
                  }
                  if (!missing(ordered_variable)) {
                    if (!is.character(ordered_variable)) {
                      stop("Argument 'ordered_variable' is not a 'character'.")
                    }
                    if (!all(ordered_variable %in% colnames(data))) {
                      stop("Some element in argument 'ordered_variable' is not recognized.")
                    } else {
                      data[, ordered_variable] <- 
                        lapply(X = ordered_variable, 
                               FUN = function(ordered_variable_i) {
                                 return(as.ordered(getElement(data, ordered_variable_i)))
                               })
                    }
                  }
                  numeric_variable <- 
                    colnames(data)[sapply(X = data, FUN = is.numeric)]
                  ordered_variable <- 
                    colnames(data)[sapply(X = data, FUN = is.ordered)]
                  if (length(ordered_variable) > 0) {
                    nlevel_ordered <- sapply(X = data[, ordered_variable],
                                             FUN = nlevels)
                    names(nlevel_ordered) <- ordered_variable
                  } else {
                    nlevel_ordered <- numeric(0)
                  }
                  if (missing(weight_variable)) {
                    weight_variable <- character(0)
                  } else {
                    if (!is.character(weight_variable)) {
                      stop("Argument 'weight_variable' is not a 'character'.")
                    }
                    if (length(weight_variable) > 1) {
                      stop("Argument 'weight_variable' can be only of length one.")
                    }
                    if (!(weight_variable %in% colnames(data))) {
                      stop("Argument 'weight_variable' is not recognized.")
                    }
                  }
                  if (missing(auxiliary_variable)) {
                    auxiliary_variable <- character(0)
                  } else {
                    if (!is.character(auxiliary_variable)) {
                      stop("Argument 'auxiliary_variable' is not a 'character'.")
                    }
                    if (!all(auxiliary_variable %in% colnames(data))) {
                      stop("Some element in argument 'auxiliary_variable' is not recognized.")
                    }
                  }
                  if (missing(group_variable)) {
                    group_variable <- character(0)
                    level_group <- "g"
                  } else {
                    if (!is.character(group_variable)) {
                      stop("Argument 'group_variable' is not a 'character'.")
                    }
                    if (length(group_variable) > 1) {
                      stop("Argument `group_variable` can be only of length one.")
                    }
                    if (!(group_variable %in% colnames(data))) {
                      stop("Argument 'group_variable' is not recognized.")
                    }
                    level_group <-
                      sort(levels(factor(
                        getElement(data, group_variable)
                      )))
                  }
                } else {
                  if (!is.matrix(sample_cov) & !is.list(sample_cov)) {
                    stop(
                      "Argument 'sample_cov' must be a 'matrix' (for single group analysis)",
                      " or a 'list' of 'matrix' (for multiple group analysis)."
                    )
                  }
                  if (is.matrix(sample_cov)) {
                    sample_cov <- list(sample_cov)
                  }
                  if (is.null(names(sample_cov))) {
                    if (length(sample_cov) == 1L) {
                      level_group <- "g"
                    } else {
                      level_group <- paste0("g", 1:length(sample_cov))
                    }
                    names(sample_cov) <- level_group
                  } else {
                    level_group <- names(sample_cov)
                  }
                  if (!missing(numeric_variable)) {
                    stop("Argument 'numeric_variable' is unnecessary under moment initialization.")
                  }
                  if (!missing(ordered_variable)) {
                    stop("Argument 'ordered_variable' is unnecessary under moment initialization.")
                  }
                  if (!missing(weight_variable)) {
                    stop("Argument 'weight_variable' is unnecessary under moment initialization.")
                  }
                  if (!missing(auxiliary_variable)) {
                    stop("Argument 'auxiliary_variable' is unnecessary under moment initialization.")
                  }
                  if (!missing(group_variable)) {
                    stop("Argument 'group_variable' is unnecessary under moment initialization.")
                  }
                  numeric_variable <- colnames(sample_cov[[1]])
                  ordered_variable <- character(0)
                  weight_variable <- character(0)
                  auxiliary_variable <- character(0)
                  group_variable <- character(0)
                  nlevel_ordered <- numeric(0)
                } 
              }
              if (any(grepl(pattern = "/|\\||@",
                            x = level_group))) {
                stop(
                  "Levels of groups cannot contain '/', '|', and '@'.",
                  "\n  Please modify the levels of groups in the specified data source."
                )
              }
              if (missing(reference_group)) {
                reference_group <- NULL
              } else {
                if (!anyNA(reference_group)) {
                  if (!is.null(reference_group)) {
                    if (length(level_group) == 1L) {
                      stop("Argument 'reference_group' is unnecessary for single group analysis.")
                    } else {
                      if (!(reference_group %in% level_group)) {
                        stop(
                          "Argument 'reference_group' is not recognized.",
                          "\n  Group names currently recognized by 'lslx' are ",
                          do.call(paste, as.list(level_group)),
                          " (possibly automatically created).",
                          "\n  Specified 'reference_group' is ",
                          reference_group,
                          "."
                        )
                      }
                    }
                  }
                }
              }
              private$model <-
                lslxModel$new(
                  model = model,
                  numeric_variable = numeric_variable,
                  ordered_variable = ordered_variable,
                  weight_variable = weight_variable,
                  auxiliary_variable = auxiliary_variable,
                  group_variable = group_variable,
                  reference_group = reference_group,
                  level_group = level_group,
                  nlevel_ordered = nlevel_ordered
                )
              private$data <-
                lslxData$new(
                  data = data,
                  sample_cov = sample_cov,
                  sample_mean = sample_mean,
                  sample_size = sample_size,
                  sample_moment_acov = sample_moment_acov,
                  numeric_variable = private$model$numeric_variable,
                  ordered_variable = private$model$ordered_variable,
                  weight_variable = private$model$weight_variable,
                  auxiliary_variable = private$model$auxiliary_variable,
                  group_variable = private$model$group_variable,
                  reference_group = private$model$reference_group,
                  level_group = private$model$level_group,
                  nlevel_ordered = private$model$nlevel_ordered,
                  name_response = private$model$name_response
                )
              private$fitting <- NULL
              if (verbose) {
                cat(
                  "An 'lslx' R6 class is initialized via",
                  ifelse(!missing(data),
                         "'data'",
                         "'sample_cov'"),
                  "argument. \n"
                )
                cat("  Response Variables:",
                    private$model$name_response,
                    "\n")
                if (length(private$model$name_factor) > 0) {
                  cat("  Latent Factors:",
                      private$model$name_factor,
                      "\n")
                }
                if (length(private$data$auxiliary) > 0) {
                  cat("  Auxiliary Variables:",
                      colnames(private$data$auxiliary[[1]]),
                      "\n")
                }
                if (length(private$model$level_group) > 1) {
                  cat("  Groups:",
                      private$model$level_group,
                      "\n")
                  if (!is.null(private$model$reference_group)) {
                    cat("  Reference Group:",
                        private$model$reference_group,
                        "\n")
                  }
                }
                if (!is.null(private$model$reference_group)) {
                  cat(
                    "NOTE:",
                    "Because",
                    private$model$reference_group,
                    "is set as reference,",
                    "coefficients in other groups actually represent increments from the reference.\n"
                  )
                }
              }
            })
