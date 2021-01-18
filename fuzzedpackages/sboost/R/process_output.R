#' @importFrom rlang .data
#' @importFrom stats sd

# --------------------------------------------------------------------------------
# PREPARES CLASSIFIER OUTPUT
process_classifier_output <- function(classifier, features, outcomes, otcm_def, call, rows_included = NULL) {

  # create classifier data frame
  clfr <- data.frame(matrix(ncol = 8, nrow = length(classifier)))
  colnames(clfr) <- c("stump", "feature", "vote", "left", "right", "split", "left_categories", "right_categories")

  # set output values
  for (i in seq_along(classifier)) {
    feature <- classifier[[i]][[1]] + 1
    orientation <- classifier[[i]][[2]]
    vote <- classifier[[i]][[3]]
    categorical <- classifier[[i]][[4]]
    split <- classifier[[i]][[5]]
    left_categories <- classifier[[i]][[6]]
    right_categories <- classifier[[i]][[7]]

    # stump
    clfr$stump[i] <- i

    # feature name
    clfr$feature[i] <- colnames(features)[[feature]]

    # vote
    clfr$vote[i] <- vote

    if (categorical == 0) {
      # orientation
      if (orientation == 1) {
        clfr$left[i] <- otcm_def["negative"]
        clfr$right[i] <- otcm_def["positive"]
      } else {
        clfr$left[i] <- otcm_def["positive"]
        clfr$right[i] <- otcm_def["negative"]
      }
      # split
      clfr$split[i] <- split
      clfr$left_categories[i] <- NA
      clfr$right_categories[i] <- NA
    }
    if (categorical == 1) {
      # orientation
      clfr$left[i] <- otcm_def["positive"]
      clfr$right[i] <- otcm_def["negative"]

      feature_levels <- levels(addNA(factor(features[[feature]])))
      # left_categories
      temp_categories <- rep(NA, length(left_categories))
      for (j in 1:length(left_categories)) {
        temp_categories[[j]] <- feature_levels[[left_categories[[j]]]]
      }
      clfr$left_categories[i] <- paste(temp_categories, collapse = "; ")
      # right_categories
      temp_categories <- rep(NA, length(right_categories))
      for (j in 1:length(right_categories)) {
        temp_categories[[j]] <- feature_levels[[right_categories[[j]]]]
      }
      clfr$right_categories[i] <- paste(temp_categories, collapse = "; ")

      # split
      clfr$split[i] <- NA
    }
  }

  # Training set information
  if (is.null(rows_included)) {
    training = data.frame(stumps = nrow(clfr),
                          features = ncol(features),
                          instances = nrow(features),
                          positive_prevalence = sum(outcomes == otcm_def["positive"]) / length(outcomes))
  } else {
    # if rows_included is not null then only a portion of the instances were used in training
    training = data.frame(stumps = nrow(clfr),
                          features = ncol(features),
                          instances = nrow(features[rows_included, ]),
                          positive_prevalence = sum(outcomes[rows_included] == otcm_def["positive"]) / length(outcomes[rows_included]))
  }

  # Assessment
  output <- list(classifier = clfr, outcomes = otcm_def, training = training, call = call)
  class(output) <- "sboost_classifier"

  return(output)
}

#' @export
print.sboost_classifier <- function(x, ...) {
  cat("SBOOST CLASSIFIER SUMMARY\n")
  cat(" ----------------------- \n")
  cat("Number of stumps trained: ", x$training$stumps, "\n", sep = "")
  cat("Number of training features: ", x$training$features, "\n\n", sep = "")

  cat("Number of training instances: ", x$training$instances, "\n", sep = "")
  cat("Positive outcome: ", as.character(x$outcomes["positive"]), "\n", sep = "")
  cat("Positive prevalence: ", x$training$positive_prevalence, "\n", sep = "")
  cat("Negative outcome: ", as.character(x$outcomes["negative"]), "\n", sep = "")
}






# --------------------------------------------------------------------------------
# PREPARE ASSESSMENT OUTPUT
process_assessment_output <- function(cumulative_statistics, feature_scores, object, call) {
  # Performance
  performance <- c(
    "true_positive" = cumulative_statistics$true_positive[nrow(cumulative_statistics)],
    "false_negative" = cumulative_statistics$false_negative[nrow(cumulative_statistics)],
    "true_negative" = cumulative_statistics$true_negative[nrow(cumulative_statistics)],
    "false_positive" = cumulative_statistics$false_positive[nrow(cumulative_statistics)],
    "prevalence" = cumulative_statistics$prevalence[nrow(cumulative_statistics)],
    "accuracy" = cumulative_statistics$accuracy[nrow(cumulative_statistics)],
    "sensitivity" = cumulative_statistics$sensitivity[nrow(cumulative_statistics)],
    "specificity" = cumulative_statistics$specificity[nrow(cumulative_statistics)],
    "ppv" = cumulative_statistics$ppv[nrow(cumulative_statistics)],
    "npv" = cumulative_statistics$npv[nrow(cumulative_statistics)],
    "f1" = cumulative_statistics$f1[nrow(cumulative_statistics)])

  # Output
  if (is.null(feature_scores)) {
    output <- list(performance = performance,
                   cumulalative_statistics = cumulative_statistics,
                   sboost_classifier = object,
                   outcomes = object$outcomes,
                   call = call)
  } else {
    output <- list(performance = performance,
                   cumulalative_statistics = cumulative_statistics,
                   feature_scores = feature_scores,
                   sboost_classifier = object,
                   outcomes = object$outcomes,
                   call = call)
  }
  class(output) <- "sboost_assessment"

  return(output)
}

#' @export
print.sboost_assessment <- function(x, ...) {
  cat("SBOOST ASSESSMENT SUMMARY\n", sep = "")
  cat(" ----------------------- \n", sep = "")
  cat("Number of classifier stumps used: ", nrow(x$cumulative_statistics), "\n", sep = "")
  cat("Number of instances assessed: ", x$performance["true_positive"] + x$performance["false_negative"] + x$performance["true_negative"] + x$performance["false_positive"], "\n\n", sep = "")

  cat("Positive outcome: ", as.character(x$outcomes["positive"]), "\n", sep = "")
  cat("Prevalence: ", x$performance["prevalence"], "\n", sep = "")
  cat("Negative outcome: ", as.character(x$outcomes["negative"]), "\n\n", sep = "")

  cat("Accuracy: ", x$performance["accuracy"], "\n", sep = "")
  cat("Sensitivity: ", x$performance["sensitivity"], "\n", sep = "")
  cat("Specificity: ", x$performance["specificity"], "\n", sep = "")
  cat("PPV: ", x$performance["ppv"], "\n", sep = "")
  cat("NPV: ", x$performance["npv"], "\n", sep = "")
  cat("F1: ", x$performance["f1"], "\n", sep = "")
}







# --------------------------------------------------------------------------------
# PREPARE VALIDATION OUTPUT
process_validation_output <- function(training_statistics, testing_statistics, classifier_list, k_fold, call) {

  last_row <- nrow(training_statistics[[1]])
  training_summary_statistics <- training_statistics[[1]]
  testing_summary_statistics <- testing_statistics[[1]]
  for (i in seq_along(training_statistics)) {
    training_summary_statistics <- rbind(training_summary_statistics, training_statistics[[i]])
    testing_summary_statistics <- rbind(testing_summary_statistics, testing_statistics[[i]])
  }

  training_summary_statistics <- dplyr::group_by(training_summary_statistics, .data$next_stump)
  training_summary_statistics <- dplyr::summarise(training_summary_statistics,
                                                  true_positive_mean = mean(.data$true_positive),
                                                  true_positive_sd = sd(.data$true_positive),
                                                  false_negative_mean = mean(.data$false_negative),
                                                  false_negative_sd = sd(.data$false_negative),
                                                  true_negative_mean = mean(.data$true_negative),
                                                  true_negative_sd = sd(.data$true_negative),
                                                  false_positive_mean = mean(.data$false_positive),
                                                  false_positive_sd = sd(.data$false_positive),
                                                  accuracy_mean = mean(.data$accuracy),
                                                  accuracy_sd = sd(.data$accuracy),
                                                  sensitivity_mean = mean(.data$sensitivity),
                                                  sensitivity_sd = sd(.data$sensitivity),
                                                  specificity_mean = mean(.data$specificity),
                                                  specificity_sd = sd(.data$specificity),
                                                  ppv_mean = mean(.data$ppv),
                                                  ppv_sd = sd(.data$ppv),
                                                  npv_mean = mean(.data$npv),
                                                  npv_sd = sd(.data$npv),
                                                  f1_mean = mean(.data$f1),
                                                  f1_sd = sd(.data$f1))
  testing_summary_statistics <- dplyr::group_by(testing_summary_statistics, .data$next_stump)
  testing_summary_statistics <- dplyr::summarise(testing_summary_statistics,
                                                 true_positive_mean = mean(.data$true_positive),
                                                 true_positive_sd = sd(.data$true_positive),
                                                 false_negative_mean = mean(.data$false_negative),
                                                 false_negative_sd = sd(.data$false_negative),
                                                 true_negative_mean = mean(.data$true_negative),
                                                 true_negative_sd = sd(.data$true_negative),
                                                 false_positive_mean = mean(.data$false_positive),
                                                 false_positive_sd = sd(.data$false_positive),
                                                 accuracy_mean = mean(.data$accuracy),
                                                 accuracy_sd = sd(.data$accuracy),
                                                 sensitivity_mean = mean(.data$sensitivity),
                                                 sensitivity_sd = sd(.data$sensitivity),
                                                 specificity_mean = mean(.data$specificity),
                                                 specificity_sd = sd(.data$specificity),
                                                 ppv_mean = mean(.data$ppv),
                                                 ppv_sd = sd(.data$ppv),
                                                 npv_mean = mean(.data$npv),
                                                 npv_sd = sd(.data$npv),
                                                 f1_mean = mean(.data$f1),
                                                 f1_sd = sd(.data$f1))

  performance <- data.frame(
    row.names = c("true_positive", "false_negative", "true_negative", "false_positive", "accuracy", "sensitivity", "specificity", "ppv", "npv", "f1"),
    training_mean = c(training_summary_statistics$true_positive_mean[last_row],
                      training_summary_statistics$false_negative_mean[last_row],
                      training_summary_statistics$true_negative_mean[last_row],
                      training_summary_statistics$false_positive_mean[last_row],
                      training_summary_statistics$accuracy_mean[last_row],
                      training_summary_statistics$sensitivity_mean[last_row],
                      training_summary_statistics$specificity_mean[last_row],
                      training_summary_statistics$ppv_mean[last_row],
                      training_summary_statistics$npv_mean[last_row],
                      training_summary_statistics$f1_mean[last_row]),
    training_sd = c(training_summary_statistics$true_positive_sd[last_row],
                    training_summary_statistics$false_negative_sd[last_row],
                    training_summary_statistics$true_negative_sd[last_row],
                    training_summary_statistics$false_positive_sd[last_row],
                    training_summary_statistics$accuracy_sd[last_row],
                    training_summary_statistics$sensitivity_sd[last_row],
                    training_summary_statistics$specificity_sd[last_row],
                    training_summary_statistics$ppv_sd[last_row],
                    training_summary_statistics$npv_sd[last_row],
                    training_summary_statistics$f1_sd[last_row]),
    testing_mean = c(testing_summary_statistics$true_positive_mean[last_row],
                     testing_summary_statistics$false_negative_mean[last_row],
                     testing_summary_statistics$true_negative_mean[last_row],
                     testing_summary_statistics$false_positive_mean[last_row],
                     testing_summary_statistics$accuracy_mean[last_row],
                     testing_summary_statistics$sensitivity_mean[last_row],
                     testing_summary_statistics$specificity_mean[last_row],
                     testing_summary_statistics$ppv_mean[last_row],
                     testing_summary_statistics$npv_mean[last_row],
                     testing_summary_statistics$f1_mean[last_row]),
    testing_sd = c(testing_summary_statistics$true_positive_sd[last_row],
                   testing_summary_statistics$false_negative_sd[last_row],
                   testing_summary_statistics$true_negative_sd[last_row],
                   testing_summary_statistics$false_positive_sd[last_row],
                   testing_summary_statistics$accuracy_sd[last_row],
                   testing_summary_statistics$sensitivity_sd[last_row],
                   testing_summary_statistics$specificity_sd[last_row],
                   testing_summary_statistics$ppv_sd[last_row],
                   testing_summary_statistics$npv_sd[last_row],
                   testing_summary_statistics$f1_sd[last_row]))

  output <- list(performance = performance,
                 training_summary_statistics = training_summary_statistics,
                 testing_summary_statistics = testing_summary_statistics,
                 training_statistics = training_statistics,
                 testing_statistics = testing_statistics,
                 classifier_list = classifier_list,
                 outcomes = classifier_list[[1]]$outcomes,
                 k_fold = k_fold,
                 call = call)
  class(output) <- "sboost_validation"
  return(output)

}

#' @export
print.sboost_validation <- function(x, ...) {
  cat("SBOOST VALIDATION SUMMARY\n", sep = "")
  cat(" ----------------------- \n", sep = "")
  cat("Number of validation sets: ", x$k_fold, "\n", sep = "")
  cat("Number of training features: ", x$classifier_list[[1]]$training$features, "\n", sep = "")
  cat("Number of stumps per training set: ", x$classifier_list[[1]]$training$stumps, "\n", sep = "")

  cat("Approximate instances per training set: ", x$classifier_list[[1]]$training$instances, "\n", sep = "")
  cat("Approximate instances per testing set: ", floor(x$classifier_list[[1]]$training$instances / (x$k_fold - 1)), "\n", sep = "")
  cat("Positive outcome: ", as.character(x$outcomes["positive"]), "\n", sep = "")
  cat("Negative outcome: ", as.character(x$outcomes["negative"]), "\n\n", sep = "")

  cat("Mean accuracy (SD)\n")
  cat(" Training: ", x$performance["accuracy", "training_mean"], " (", x$performance["accuracy", "training_sd"], ")\n", sep = "")
  cat(" Testing: ", x$performance["accuracy", "testing_mean"], " (", x$performance["accuracy", "testing_sd"], ")\n", sep = "")
  cat("Mean F1 (SD)\n")
  cat(" Training: ", x$performance["f1", "training_mean"], " (", x$performance["f1", "training_sd"], ")\n", sep = "")
  cat(" Testing: ", x$performance["f1", "testing_mean"], " (", x$performance["f1", "testing_sd"], ")\n", sep = "")
}






