
## \code{$new()} initializes a new \code{lslxData} object. ##
lslxData$set("public",
             "initialize",
             function(data,
                      sample_cov,
                      sample_mean,
                      sample_size,
                      sample_moment_acov,
                      numeric_variable,
                      ordered_variable,
                      weight_variable,
                      auxiliary_variable,
                      group_variable,
                      reference_group,
                      level_group,
                      nlevel_ordered,
                      name_response) {
               if (!missing(data)) {
                 index <- 1:nrow(data)
                 row.names(data) <- index
                 if (length(group_variable) == 0) {
                   response <-
                     list(data[, name_response, drop = FALSE])
                   names(response) <- level_group
                   if (length(weight_variable) == 0) {
                     weight <- 
                       list(data.frame(weight = rep(1, nrow(data)),
                                       row.names = index))
                   } else {
                     weight <- 
                       list(data[, weight_variable, drop = FALSE])
                   }
                   names(weight) <- level_group
                   if (length(auxiliary_variable) == 0) {
                     auxiliary <- list()
                   } else {
                     auxiliary <- 
                       list(data[, auxiliary_variable, drop = FALSE])
                   }
                 } else {
                   data <-
                     data[order(as.character(getElement(data, group_variable))), 
                          , 
                          drop = FALSE]
                   data[, group_variable] <-
                     as.character(getElement(data, group_variable))
                   response <-
                     split(data[, name_response, drop = FALSE],
                           getElement(data, group_variable))
                   if (length(weight_variable) == 0) {
                     weight <- 
                       split(data.frame(weight = rep(1, nrow(data)),
                                        row.names = row.names(data)), 
                             getElement(data, group_variable))
                   } else {
                     weight <-
                       split(data[, weight_variable, drop = FALSE],
                             getElement(data, group_variable))
                   }
                   if (length(auxiliary_variable) == 0) {
                     auxiliary <- list()
                   } else {
                     auxiliary <-
                       split(data[, auxiliary_variable, drop = FALSE],
                             getElement(data, group_variable))
                   }
                 }
                 if (!all(sapply(X = response, 
                                 FUN = function(response_i) {
                                   sapply(X = response_i,
                                          FUN = function(response_ij) {
                                            return(is.numeric(response_ij) |
                                                   is.ordered(response_ij))
                                          })
                                 }))) {
                   stop("Response variables cannot contain non-numeric/ordered variables.")
                 }
                 if (length(auxiliary) > 0) {
                   if (!all(sapply(X = auxiliary, 
                                   FUN = function(auxiliary_i) {
                                     sapply(X = auxiliary_i,
                                            FUN = function(auxiliary_ij) {
                                              return(is.numeric(auxiliary_ij))
                                            })
                                   }))) {
                     stop("Auxiliary variable(s) cannot contain non-numeric variables.")
                   }
                 }
                 if (any(sapply(
                   X = weight,
                   FUN = function(weight_i) {
                     return(any(weight_i < 0))
                   }
                 ))) {
                   stop("Weight variable cannot contain negative value.")
                 }
                 self$response <- response
                 self$weight <- weight
                 self$auxiliary <- auxiliary
                 self$pattern <-
                   lapply(
                     X = self$response,
                     FUN = function(response_i) {
                       return(!is.na(response_i))
                     }
                   )
                 self$sample_cov <- list()
                 self$sample_mean <- list()
                 self$sample_size <- list()
                 self$sample_moment_acov <- list()
               } else {
                 sample_cov <-
                   lapply(
                     X = sample_cov,
                     FUN = function(sample_cov_i) {
                       sample_cov_i <-
                         sample_cov_i[name_response,
                                      name_response,
                                      drop = FALSE]
                       return(sample_cov_i)
                     }
                   )
                 names(sample_cov) <- level_group
                 if (missing(sample_mean)) {
                   sample_mean <-
                     lapply(
                       X = level_group,
                       FUN = function(i) {
                         sample_mean_i <- rep(0, ncol(sample_cov[[1]]))
                         names(sample_mean_i) <-
                           colnames(sample_cov[[1]])
                         return(sample_mean_i)
                       }
                     )
                   names(sample_mean) <- level_group
                 } else {
                   if (!is.numeric(sample_mean) & !is.list(sample_mean)) {
                     stop("Argument 'sample_mean' must be 'numeric' or 'list' of 'numeric'.")
                   }
                   if (is.numeric(sample_mean)) {
                     sample_mean <- list(sample_mean)
                   }
                   if (!all(name_response %in% names(sample_mean[[1]]))) {
                     stop(
                       "Some response variable in 'model' cannot be found in 'sample_mean'.",
                       "\n  Response variables specified by 'model' are ",
                       do.call(paste, as.list(name_response)),
                       ".",
                       "\n  Column names of 'sample_mean' are ",
                       do.call(paste, as.list(colnames(
                         sample_mean[[1]]
                       ))),
                       "."
                     )
                   } else {
                     sample_mean <-
                       lapply(
                         X = sample_mean,
                         FUN = function(sample_mean_i) {
                           sample_mean_i <-
                             sample_mean_i[name_response]
                           return(sample_mean_i)
                         }
                       )
                   }
                   if (length(sample_mean) != length(sample_cov)) {
                     stop(
                       "The length of argument 'sample_mean' doesn't match the length of 'sample_cov'.",
                       "\n  The length of 'sample_mean' is ",
                       length(sample_mean),
                       ".",
                       "\n  The length of 'sample_cov' is ",
                       length(sample_cov),
                       "."
                     )
                   }
                   if (is.null(names(sample_mean))) {
                     names(sample_mean) <- level_group
                   } else {
                     if (!all(level_group %in% names(sample_mean))) {
                       stop(
                         "Argument 'sample_mean' contains unrecognized group name(s).",
                         "\n  Group name(s) currently recognized by 'lslx' is ",
                         do.call(paste, as.list(level_group)),
                         " (possibly automatically created).",
                         "\n  Group name(s) specified in 'sample_mean' is ",
                         do.call(paste, as.list(names(
                           sample_mean
                         ))),
                         "."
                       )
                     } else{
                       sample_mean <- sample_mean[level_group]
                     }
                   }
                 }
                 if (missing(sample_size)) {
                   stop("Argument 'sample_size' cannot be empty if 'sample_cov' is used.")
                 } else {
                   if (!is.numeric(sample_size) & !is.list(sample_size)) {
                     stop("Argument 'sample_size' must be a 'numeric' or a 'list' of 'numeric'.")
                   }
                   if (is.numeric(sample_size)) {
                     sample_size <- as.list(sample_size)
                   }
                   if (length(sample_size) != length(sample_cov)) {
                     stop(
                       "The length of argument 'sample_size' doesn't match the length of 'sample_cov'.",
                       "\n  The length of 'sample_size' is ",
                       length(sample_size),
                       ".",
                       "\n  The length of 'sample_cov' is ",
                       length(sample_cov),
                       "."
                     )
                   }
                   if (is.null(names(sample_size))) {
                     names(sample_size) <- level_group
                   } else {
                     if (!all(level_group %in% names(sample_size))) {
                       stop(
                         "Argument 'sample_size' contains unrecognized group name(s).",
                         "\n  Group name(s) currently recognized by 'lslx' is ",
                         do.call(paste, as.list(level_group)),
                         " (possibly automatically created).",
                         "\n  Group name(s) specified in 'sample_size' is ",
                         do.call(paste, as.list(names(
                           sample_size
                         ))),
                         "."
                       )
                     } else{
                       sample_size <- sample_size[level_group]
                     }
                   }
                 }
                 
                 if (!missing(sample_moment_acov)) {
                   if (!is.matrix(sample_moment_acov) & !is.list(sample_moment_acov)) {
                     stop(
                       "Argument 'sample_moment_acov' must be a 'matrix' (for single group analysis)",
                       " or a 'list' of 'matrix' (for multiple group analysis)."
                     )
                   }
                   if (is.matrix(sample_moment_acov)) {
                     sample_moment_acov <- list(sample_moment_acov)
                   }
                   
                   name_response2 <- outer(
                     name_response,
                     name_response,
                     FUN = function(name_response_i, name_response_j) {
                       return(paste(name_response_i, name_response_j, sep = "*"))
                     }
                   )
                   name_response2 <-
                     name_response2[lower.tri(name_response2, diag = TRUE)]
                   if (!all(c(name_response, name_response2) %in% 
                            colnames(sample_moment_acov[[1]]))) {
                     stop(
                       "Some response variables (products) in 'model' cannot be found in 'sample_moment_acov'.",
                       "\n  Response variables (products) specified by 'model' are ",
                       do.call(paste, as.list(c(name_response, name_response2))),
                       ".",
                       "\n  Column names of 'sample_moment_acov' are ",
                       do.call(paste, as.list(colnames(
                         sample_moment_acov[[1]]
                       ))),
                       "."
                     )
                   } else {
                     sample_moment_acov <-
                       lapply(
                         X = sample_moment_acov,
                         FUN = function(sample_moment_acov_i) {
                           sample_moment_acov_i <-
                             sample_moment_acov_i[c(name_response, name_response2),
                                                  c(name_response, name_response2),
                                                  drop = FALSE]
                           return(sample_moment_acov_i)
                         }
                       )
                   }
                   if (length(sample_moment_acov) != length(sample_cov)) {
                     stop(
                       "The length of argument 'sample_moment_acov' doesn't match the length of 'sample_cov'.",
                       "\n  The length of 'sample_moment_acov' is ",
                       length(sample_moment_acov),
                       ".",
                       "\n  The length of 'sample_cov' is ",
                       length(sample_cov),
                       "."
                     )
                   }
                   if (is.null(names(sample_moment_acov))) {
                     names(sample_moment_acov) <- level_group
                   } else {
                     if (!all(level_group %in% names(sample_moment_acov))) {
                       stop(
                         "Argument 'sample_moment_acov' contains unrecognized group name(s).",
                         "\n  Group name(s) currently recognized by 'lslx' is ",
                         do.call(paste, as.list(level_group)),
                         " (possibly automatically created).",
                         "\n  Group name(s) specified in 'sample_moment_acov' is ",
                         do.call(paste, as.list(names(
                           sample_moment_acov
                         ))),
                         "."
                       )
                     } else{
                       sample_moment_acov <- sample_moment_acov[level_group]
                     }
                   }
                 } else {
                   sample_moment_acov <- list()
                 }
                 self$response <- list()
                 self$weight <- list()
                 self$auxiliary <- list()
                 self$pattern <- list()
                 self$sample_cov <- sample_cov
                 self$sample_mean <- sample_mean
                 self$sample_size <- sample_size
                 self$sample_moment_acov <- sample_moment_acov
               }
               self$numeric_variable <- numeric_variable
               self$ordered_variable <- ordered_variable
               self$weight_variable <- weight_variable
               self$auxiliary_variable <- auxiliary_variable
               self$group_variable <- group_variable
               self$reference_group <- reference_group
               self$level_group <- level_group
               self$nlevel_ordered <- nlevel_ordered
               self$name_response <- name_response
             })

