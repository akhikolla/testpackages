## \code{$plot_numerical_condition()} plots the values of selected numerical conditions. ##
lslx$set("public",
         "plot_numerical_condition",
         function(condition,
                  x_scale = "default",
                  x_reverse = "default",
                  mode = "default") {
           if (private$fitting$control$regularizer) {
             if (length(private$fitting$control$lambda_grid) <= 1) {
               stop(
                 "When 'length(lambda_grid) == 1', 'plot_numerical_condition()' method is not available."
               )
             }
           } else {
             
           }
           if (private$fitting$control$searcher) {
             if (length(private$fitting$control$step_grid) <= 1) {
               stop(
                 "When 'length(step_grid) == 1', 'plot_numerical_condition()' method is not available."
               )
             }
           } else {
             
           }
           if (private$fitting$control$penalty_method == "none") {
             stop(
               "When 'penalty_method' is 'none', 'plot_numerical_condition()' method is not available."
             )
           }
           if (private$fitting$control$double_regularizer) {
             stop(
               "When two regularizers are implemented, 'plot_numerical_condition()' method is not available."
             )
           }
           if (missing(condition)) {
             condition <-
               c("n_iter_out",
                 "objective_gradient_abs_max",
                 "objective_hessian_convexity")
           } else {
             if (condition == "all") {
               condition <-
                 setdiff(
                   names(private$fitting$fitted_result$numerical_condition[[1]]),
                   c("lambda_1st", "lambda_2nd", "delta_1st", "delta_2nd", "step")
                 )
             }
             if (any(!(
               condition %in% names(private$fitting$fitted_result$numerical_condition[[1]])
             ))) {
               stop("Argument 'condition' contains unrecognized numerical condition.")
             }
           }
           if (x_scale == "default") {
             x_scale <- "identity"
           }
           if (x_reverse == "default") {
             if (private$fitting$control$penalty_method == "forward") {
               x_reverse <- TRUE
             } else {
               x_reverse <- FALSE
             }
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$numerical_condition)[condition,
                                                                                      ,
                                                                                      drop = FALSE])
           condition <-
             ifelse(
               condition == "objective_gradient_abs_max",
               "gradient",
               ifelse(
                 condition == "objective_hessian_convexity",
                 "hessian",
                 condition
               )
             )
           condition <-
             gsub(pattern = "_",
                  replacement = " ",
                  x = condition)
           
           df_for_plot$condition <- condition
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "condition",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|\\,|\\(|\\)")
           if (private$fitting$control$regularizer) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$lambda <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[3]
                 }
               ))
             df_for_plot$delta <-
               round(as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[8]
                 }
               )), 3)
           } else {
             
           }
           if (private$fitting$control$searcher) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$step <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[2]
                 }
               ))
           } else {
             
           }
           if (mode == "plot") {
             p <- ggplot2::ggplot(df_for_plot)
             if (private$fitting$control$regularizer) {
               p <- p + ggplot2::geom_line(ggplot2::aes(x = lambda, y = value))
             } else {
               
             }
             if (private$fitting$control$searcher) {
               p <- p + ggplot2::geom_step(ggplot2::aes(x = step, y = value))
             } else {
               
             }
             if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
               p <- p + ggplot2::facet_grid(
                 condition ~ delta,
                 scales = "free_y",
                 labeller = ggplot2::labeller(
                   delta = ggplot2::label_both,
                   condition = ggplot2::label_value
                 )
               )
             } else {
               p <- p + ggplot2::facet_grid(
                 condition ~ .,
                 scales = "free_y",
                 labeller = ggplot2::labeller(
                   delta = ggplot2::label_both,
                   condition = ggplot2::label_value
                 )
               )
             }
             if (x_reverse) {
               p <- p + ggplot2::scale_x_reverse()
             }
             p + ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             ) +
               ggplot2::coord_trans(x = x_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Numerical Conditions across Penalty Levels"),
                 x = ifelse(private$fitting$control$regularizer,
                            "lambda", "step"),
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_information_criterion()} shows how the values of information criteria vary with penalty levels. ##
lslx$set("public",
         "plot_information_criterion",
         function(criterion,
                  x_scale = "default",
                  x_reverse = "default",
                  mode = "default") {
           if (private$fitting$control$regularizer) {
             if (length(private$fitting$control$lambda_grid) <= 1) {
               stop(
                 "When 'length(lambda_grid) == 1', 'plot_information_criterion()' method is not available."
               )
             }
           } else {
             
           }
           if (private$fitting$control$searcher) {
             if (length(private$fitting$control$step_grid) <= 1) {
               stop(
                 "When 'length(step_grid) == 1', 'plot_information_criterion()' method is not available."
               )
             }
           } else {
             
           }
           if (private$fitting$control$penalty_method == "none") {
             stop(
               "When 'penalty_method' is 'none', 'plot_information_criterion()' method is not available."
             )
           }
           if (private$fitting$control$double_regularizer) {
             stop(
               "When two regularizers are implemented, 'plot_information_criterion()' method is not available."
             )
           }
           if (missing(criterion)) {
             criterion <- c("aic", "aic3", "caic", "bic", "abic", "hbic")
           } else {
             if (any(!(
               criterion %in% names(private$fitting$fitted_result$information_criterion[[1]])
             ))) {
               stop("Argument `criterion` contains unrecognized information criterion.")
             }
           }
           if (x_scale == "default") {
             x_scale <- "identity"
           }
           if (x_reverse == "default") {
             if (private$fitting$control$penalty_method == "forward") {
               x_reverse <- TRUE
             } else {
               x_reverse <- FALSE
             }
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(
               cbind,
               private$fitting$fitted_result$information_criterion
             )[criterion, , drop = FALSE])
           df_for_plot$criterion <- criterion
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "criterion",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|\\,|\\(|\\)")
           if (private$fitting$control$regularizer) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$lambda <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[3]
                 }
               ))
             df_for_plot$delta <-
               round(as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[8]
                 }
               )), 3)
           } else {
             
           }
           
           if (private$fitting$control$searcher) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$step <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[2]
                 }
               ))
           } else {
             
           }
           if (mode == "plot") {
             p <- ggplot2::ggplot(df_for_plot)
             if (private$fitting$control$regularizer) {
               p <-
                 p + ggplot2::geom_line(mapping = ggplot2::aes(
                   x = lambda,
                   y = value,
                   colour = criterion
                 ))
             } else {
               
             }
             if (private$fitting$control$searcher) {
               p <-
                 p + ggplot2::geom_step(mapping = ggplot2::aes(
                   x = step,
                   y = value,
                   colour = criterion
                 ))
             } else {
               
             }
             if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
               p <-
                 p + ggplot2::facet_grid(. ~ delta, labeller = ggplot2::label_both)
             }
             if (x_reverse) {
               p <- p + ggplot2::scale_x_reverse()
             }
             p + ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             ) +
               ggplot2::coord_trans(x = x_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0("Values of Information Criteria across Penalty Levels"),
                 x = ifelse(private$fitting$control$regularizer,
                            "lambda", "step"),
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_fit_index()} shows how the values of fit indices vary with penalty levels. ##
lslx$set("public",
         "plot_fit_index",
         function(index,
                  x_scale = "default",
                  x_reverse = "default",
                  mode = "default") {
           if (private$fitting$control$regularizer) {
             if (length(private$fitting$control$lambda_grid) <= 1) {
               stop("When 'length(lambda_grid) == 1', 'plot_fit_index()' method is not available.")
             }
           } else {
             
           }
           if (private$fitting$control$searcher) {
             if (length(private$fitting$control$step_grid) <= 1) {
               stop("When 'length(step_grid) == 1', 'plot_fit_index()' method is not available.")
             }
           } else {
             
           }
           if (private$fitting$control$penalty_method == "none") {
             stop("When 'penalty_method' is 'none', 'plot_fit_index()' method is not available.")
           }
           if (private$fitting$control$double_regularizer) {
             stop(
               "When two regularizers are implemented, 'plot_fit_index()' method is not available."
             )
           }
           if (missing(index)) {
             index <-
               names(private$fitting$fitted_result$fit_index[[1]])
           } else {
             if (any(!(
               index %in% names(private$fitting$fitted_result$fit_index[[1]])
             ))) {
               stop("Argument `index` contains unrecognized fit index.")
             }
           }
           if (x_scale == "default") {
             x_scale <- "identity"
           }
           if (x_reverse == "default") {
             if (private$fitting$control$penalty_method == "forward") {
               x_reverse <- TRUE
             } else {
               x_reverse <- FALSE
             }
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$fit_index)[index,
                                                                            ,
                                                                            drop = FALSE])
           df_for_plot$index <- index
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "index",
               v.names = "value",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-ncol(df_for_plot)],
               times = colnames(df_for_plot)[-ncol(df_for_plot)],
               direction = "long"
             )
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|\\,|\\(|\\)")
           if (private$fitting$control$regularizer) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$lambda <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[3]
                 }
               ))
             df_for_plot$delta <-
               round(as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[8]
                 }
               )), 3)
           } else {
             
           }
           
           if (private$fitting$control$searcher) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$step <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[2]
                 }
               ))
           } else {
             
           }
           
           if (mode == "plot") {
             p <- ggplot2::ggplot(df_for_plot)
             if (private$fitting$control$regularizer) {
               p <-
                 p + ggplot2::geom_line(mapping = ggplot2::aes(
                   x = lambda,
                   y = value,
                   colour = index
                 ))
             } else {
               
             }
             if (private$fitting$control$searcher) {
               p <-
                 p + ggplot2::geom_step(mapping = ggplot2::aes(
                   x = step,
                   y = value,
                   colour = index
                 ))
             } else {
               
             }
             if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
               p <-
                 p + ggplot2::facet_grid(. ~ delta, labeller = ggplot2::label_both)
             }
             if (x_reverse) {
               p <- p + ggplot2::scale_x_reverse()
             }
             p + ggplot2::theme(
               panel.grid.minor = ggplot2::element_line(size = .1),
               panel.grid.major = ggplot2::element_line(size = .2)
             ) +
               ggplot2::coord_trans(x = x_scale, y = "identity") +
               ggplot2::ylim(0, 1) +
               ggplot2::labs(
                 title = paste0("Values of Goodness-of-Fit across Penalty Levels"),
                 x = ifelse(private$fitting$control$regularizer,
                            "lambda", "step"),
                 y = "value"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })

## \code{$plot_coefficient()} visualizes the solution paths of coefficients. ##
lslx$set("public",
         "plot_coefficient",
         function(block,
                  left,
                  right,
                  both,
                  x_scale = "default",
                  x_reverse = "default",
                  mode = "default") {
           if (private$fitting$control$regularizer) {
             if (length(private$fitting$control$lambda_grid) <= 1) {
               stop("When 'length(lambda_grid) == 1', 'plot_coefficient()' method is not available.")
             }
           } else {
             
           }
           
           if (private$fitting$control$searcher) {
             if (length(private$fitting$control$step_grid) <= 1) {
               stop("When 'length(step_grid) == 1', 'plot_coefficient()' method is not available.")
             }
           } else {
             
           }
           if (private$fitting$control$penalty_method == "none") {
             stop("When 'penalty_method' is 'none', 'plot_coefficient()' method is not available.")
           }
           if (private$fitting$control$double_regularizer) {
             stop(
               "When two regularizers are implemented, 'plot_coefficient()' method is not available."
             )
           }
           if (x_scale == "default") {
             x_scale <- "identity"
           }
           if (x_reverse == "default") {
             if (private$fitting$control$penalty_method == "forward") {
               x_reverse <- TRUE
             } else {
               x_reverse <- FALSE
             }
           }
           if (mode == "default") {
             mode <- "plot"
           } else {
             if (!(mode %in% c("plot", "return"))) {
               stop("Argument 'mode' can be only either 'default', 'plot' or 'return'. ")
             }
           }
           if (missing(block)) {
             block <- unique(private$model$specification$block)
           }
           if (missing(left)) {
             left <- c(private$model$name_eta, "1")
           }
           if (missing(right)) {
             right <- c(private$model$name_eta, "1")
           }
           if (missing(both)) {
             both <- c(private$model$name_eta, "1")
           }
           df_for_plot <-
             as.data.frame(do.call(cbind,
                                   private$fitting$fitted_result$coefficient))
           df_for_plot$name <- rownames(df_for_plot)
           df_for_plot$type <- private$model$specification$type
           df_for_plot <-
             df_for_plot[(private$model$specification$block %in% block) &
                           (private$model$specification$left %in% left) &
                           (private$model$specification$right %in% right) &
                           (private$model$specification$left %in% both) &
                           (private$model$specification$right %in% both), , drop = FALSE]
           if (nrow(df_for_plot) == 0) {
             stop("No such type of coefficient in the specified model.")
           }
           df_for_plot <-
             reshape(
               data = df_for_plot,
               idvar = "name",
               v.names = "estimate",
               timevar = "penalty_level",
               varying = colnames(df_for_plot)[-c(ncol(df_for_plot) - 1, ncol(df_for_plot))],
               times = colnames(df_for_plot)[-c(ncol(df_for_plot) - 1, ncol(df_for_plot))],
               direction = "long"
             )
           name_split <- strsplit(x = df_for_plot$name,
                                  split = "/")
           df_for_plot$name <- NULL
           df_for_plot$relation <-
             sapply(
               X = name_split,
               FUN = function(x) {
                 x[1]
               }
             )
           df_for_plot$group <-
             sapply(
               X = name_split,
               FUN = function(x) {
                 x[2]
               }
             )
           penalty_level_split <-
             strsplit(x = df_for_plot$penalty_level,
                      split = "=|\\,|\\(|\\)")
           if (private$fitting$control$regularizer) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$lambda <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[3]
                 }
               ))
             df_for_plot$delta <-
               round(as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[8]
                 }
               )), 3)
           } else {
             
           }
           
           if (private$fitting$control$searcher) {
             df_for_plot$penalty_level <- NULL
             df_for_plot$step <-
               as.numeric(sapply(
                 X = penalty_level_split,
                 FUN = function(x) {
                   x[2]
                 }
               ))
           } else {
             
           }
           if (mode == "plot") {
             p <- ggplot2::ggplot(df_for_plot)
             if (private$fitting$control$regularizer) {
               p <-
                 p + ggplot2::geom_line(
                   mapping = ggplot2::aes(
                     x = lambda,
                     y = estimate,
                     colour = relation,
                     linetype = type
                   ),
                   show.legend = c(colour = T, linetype = F)
                 )
             } else {
               
             }
             if (private$fitting$control$searcher) {
               p <-
                 p + ggplot2::geom_step(
                   mapping = ggplot2::aes(
                     x = step,
                     y = estimate,
                     colour = relation,
                     linetype = type
                   ),
                   show.legend = c(colour = T, linetype = F)
                 )
             } else {
               
             }
             if (private$fitting$reduced_model$n_group > 1) {
               if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
                 p <-
                   p + ggplot2::facet_grid(group ~ delta, labeller = ggplot2::label_both)
               } else {
                 p <-
                   p + ggplot2::facet_grid(group ~ ., labeller = ggplot2::label_both)
               }
             } else {
               if (private$fitting$control$penalty_method %in% c("elastic_net", "mcp")) {
                 p <-
                   p + ggplot2::facet_grid(. ~ delta, labeller = ggplot2::label_both)
               } else {
                 
               }
             }
             if (x_reverse) {
               p <- p + ggplot2::scale_x_reverse()
             }
             p + ggplot2::scale_linetype_manual(values = c(
               "fixed" = "solid",
               "free" = "solid",
               "pen" = "twodash"
             )) +
               ggplot2::theme(
                 panel.grid.minor = ggplot2::element_line(size = .1),
                 panel.grid.major = ggplot2::element_line(size = .2)
               ) +
               ggplot2::coord_trans(x = x_scale, y = "identity") +
               ggplot2::labs(
                 title = paste0(
                   "Solution Paths of Coefficients in Block ",
                   do.call(paste, as.list(block))
                 ),
                 x = ifelse(private$fitting$control$regularizer,
                            "lambda", "step"),
                 y = "coefficient estimate"
               ) +
               ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
           } else {
             return(df_for_plot)
           }
         })
