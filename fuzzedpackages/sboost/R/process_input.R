# --------------------------------------------------------------------------------
# TESTS AND PREPARES FEATURES
process_feature_input <- function(features) {

  if (!is.data.frame(features)) stop("Features must be data frame.")

  # Features which are logical, character, or factors must be changed to numeric
  #   to function with the C++ code.
  for (i in seq_along(features)) {
    if (is.logical(features[[i]]) || is.character(features[[i]])) {
      features[[i]] <- factor(features[[i]])
    }

    if (is.factor(features[[i]])) {
      features[[i]] <- as.numeric(features[[i]])
    } else if (!is.numeric(features[[i]])) {
      stop(paste("Unknown data type in column ", i))
    }
  }

  return(data.matrix(features))
}



# --------------------------------------------------------------------------------
# TESTS AND PREPARES OUTCOMES
process_outcome_input <- function(outcomes, features, otcm_def) {

  if (!is.vector(outcomes)) stop("Outcomes must be data frame or vector.")
  if (length(outcomes) != nrow(features)) stop("All training examples must have an outcome.")
  if (length(unique(outcomes)) > 2) stop("Only two distinct outcomes may be assessed.")

  # Outcomes must be either 1 or -1 to function with the c++ code.
  for (i in seq_along(outcomes)) {
    if (outcomes[[i]] == otcm_def["positive"]) {
      outcomes[[i]] <- 1
    } else {
      outcomes[[i]] <- -1
    }
  }

  return(as.numeric(outcomes))
}



# --------------------------------------------------------------------------------
# TESTS AND PREPARES CLASSIFIER INPUT
process_classifier_input <- function(classifier, features) {

  if (!("sboost_classifier" %in% class(classifier))) stop("Classifier must be an output from sboost.")

  new_classifier = list()

  # Classifier must be a list of numeric vectors (one for each stump)to function
  #   with C++ code. The numeric values must comport with the numeric falues
  #   generated in process_feature_input.
  # Vectors are organized as follows:
  #   1 - feature - feature index (C++ starts at zero so one is subtracted)
  #   2 - orientation - 1 is positive on the right, -1 is positive on the left
  #   3 - vote - unchanged
  #   4 - categorical - 1 if feature is categorical, 0 if numeric
  #   5 - split position IF numeric ELSE the first category in left_categories
  #   ... - the rest of the categories in left_categories
  for (i in 1:nrow(classifier$classifier)) {
    feature <- classifier$classifier$feature[i]
    vote <- classifier$classifier$vote[i]
    orientation <- classifier$classifier$left[i]
    if (is.na(classifier$classifier$split[i])) {
      categorical <- 1
      split <- 0
      left_categories <- strsplit(classifier$classifier$left_categories[i], "; ")[[1]]
      right_categories <- strsplit(classifier$classifier$right_categories[i], "; ")[[1]]
    } else {
      categorical <- 0
      split <- classifier$classifier$split[i]
      left_categories <- NULL
      right_categories <- NULL
    }

    # Change feature
    feature <- match(feature, colnames(features))[[1]] - 1

    # Change direction
    if (categorical == 1) {
      if (orientation[[1]] == classifier$outcomes["positive"]) {
        orientation <- 1
      } else {
        orientation <- -1
      }
    } else {
      if (orientation[[1]] == classifier$outcomes["negative"]) {
        orientation <- 1
      } else {
        orientation <- -1
      }
    }

    feature_levels <- levels(addNA(factor(features[[feature + 1]])))
    # Change left_categories
    if (categorical == 1) {
      for (j in seq_along(left_categories)) {
        if (!is.na(match(left_categories[[j]], feature_levels, nomatch = NA, incomparables = NA))) {
          left_categories[[j]] <- match(left_categories[[j]], feature_levels)
        }
      }
    }

    # Change right_categories
    if (categorical == 1) {
      for (j in seq_along(right_categories)) {
        if (!is.na(match(right_categories[[j]], feature_levels, nomatch = NA, incomparables = NA))) {
          right_categories[[j]] <- match(right_categories[[j]], feature_levels)
        }
      }
    }

    new_classifier[[i]] <- list(
      as.numeric(feature),
      as.numeric(orientation),
      as.numeric(vote),
      as.numeric(categorical),
      as.numeric(split),
      suppressWarnings(as.numeric(left_categories)),
      suppressWarnings(as.numeric(right_categories))
    )

  }

  return(new_classifier)
}



# --------------------------------------------------------------------------------
# TESTS POSITIVE SPECIFICATION
# Returns defined outcome possibilities
check_positive_value <- function(outcomes, positive) {
  otcm_p <- sort(unique(outcomes))
  if (length(otcm_p) < 2) stop("There must be two distinct outcomes to use sboost.")
  if (!is.null(positive)) {
    if (!(positive %in% otcm_p)) {
      warning("'positive' variable does not match one of the outcomes. The positive value will be the first outcome in alphabetical order.")
      positive <- NULL
    } else if (positive == otcm_p[[1]]) {
      return(c(positive = otcm_p[[1]], negative = otcm_p[[2]]))
    } else {
      return(c(positive = otcm_p[[2]], negative = otcm_p[[1]]))
    }
  } else {
    return(c(positive = otcm_p[[1]], negative = otcm_p[[2]]))
  }
}



# --------------------------------------------------------------------------------
# FIND CATEGORICAL VECTOR
find_categorical <- function(features) {
  categorical <- rep(0, ncol(features))

  # Produces a vector the same length as the number of features for C++ code
  #   1 if categorical
  #   0 if numeric
  for (i in seq_along(features)) {
    if (is.logical(features[[i]]) || is.character(features[[i]]) || is.factor(features[[i]])) {
      categorical[[i]] <- length(unique(features[[i]][!is.na(features[[i]])]))
    }
  }

  return(categorical)
}


