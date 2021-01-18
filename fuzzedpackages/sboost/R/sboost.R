#' @useDynLib sboost, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
.onUnload <- function (libpath) {
  library.dynam.unload("sboost", libpath)
}

#' sboost Learning Algorithm
#'
#' A machine learning algorithm using AdaBoost on decision stumps.
#'
#' Factors and characters are treated as categorical features.
#' Missing values are supported.
#'
#' See \url{https://jadonwagstaff.github.io/sboost.html} for a description
#' of the algorithm.
#'
#' For original paper describing AdaBoost see:
#'
#' Freund, Y., Schapire, R.E.: A decision-theoretic generalization of on-line
#' learning and an application to boosting. Journal of Computer and System Sciences
#' 55(1), 119-139 (1997)
#' @param features feature set data.frame.
#' @param outcomes outcomes corresponding to the features.
#' @param iterations number of boosts.
#' @param positive the positive outcome to test for; if NULL, the first outcome in
#'                 alphabetical (or numerical) order will be chosen.
#' @param verbose If true, progress bar will be displayed in console.
#' @return An \emph{sboost_classifier} S3 object containing:
#' \describe{
#'   \item{\emph{classifier}}{\emph{stump} - the index of the decision stump}
#'     \item{}{\emph{feature} - name of the column that this stump splits on.}
#'     \item{}{\emph{vote} - the weight that this stump has on the final classifier.}
#'     \item{}{\emph{orientation} - shows how outcomes are split. If \emph{feature} is numeric
#'            shows split orientation, if \emph{feature} value is less than \emph{split} then vote
#'            is cast in favor of left side outcome, otherwise the vote is cast for the
#'            right side outcome. If \emph{feature} is categorical, vote is
#'            cast for the left side outcome if \emph{feature} value is found in
#'            \emph{left_categories}, otherwise vote is cast for right side outcome.}
#'     \item{}{\emph{split} - if \emph{feature} is numeric, the value where the decision stump
#'             splits the outcomes; otherwise, NA.}
#'     \item{}{\emph{left_categories} - if \emph{feature} is categorical, shows the \emph{feature}
#'             values that sway the vote to the left side outcome on the \emph{orientation} split;
#'             otherwise, NA.}
#'   \item{\emph{outcomes}}{Shows which outcome was considered as positive and which negative.}
#'   \item{\emph{training}}{\emph{stumps} - how many decision stumps were trained.}
#'     \item{}{\emph{features} - how many features the training set contained.}
#'     \item{}{\emph{instances} - how many instances or rows the training set contained.}
#'     \item{}{\emph{positive_prevalence} - what fraction of the training instances were positive.}
#'   \item{\emph{call}}{Shows the parameters that were used to build the classifier.}
#' }
#' @seealso
#'   \code{\link{predict.sboost_classifier}} - to get predictions from the classifier.
#'
#'   \code{\link{assess}} - to evaluate the performance of the classifier.
#'
#'   \code{\link{validate}} - to perform cross validation for the classifier training.
#' @examples
#' # malware
#' malware_classifier <- sboost(malware[-1], malware[1], iterations = 5, positive = 1)
#' malware_classifier
#' malware_classifier$classifier
#'
#' # mushrooms
#' mushroom_classifier <- sboost(mushrooms[-1], mushrooms[1], iterations = 5, positive = "p")
#' mushroom_classifier
#' mushroom_classifier$classifier
#' @export
sboost <- function(features, outcomes, iterations = 1, positive = NULL, verbose = FALSE) {

  # PREPARE INPUT
  # --------------------------------------------------------------------------------
  if (is.data.frame(outcomes)) outcomes <- as.vector(outcomes[[1]])
  processed_features <- process_feature_input(features)
  categorical <- find_categorical(features)
  otcm_def <- check_positive_value(outcomes, positive)
  processed_outcomes <- process_outcome_input(outcomes, features, otcm_def)


  # DEVELOP CLASSIFIER
  # --------------------------------------------------------------------------------
  classifier <- make_classifier(processed_features, processed_outcomes, categorical, iterations, verbose)
  classifier <- process_classifier_output(classifier, features, outcomes, otcm_def, match.call())

  return(classifier)
}



# make_classifier takes an unordered set of features and outcomes,
#        orders them, and calls the appropriate functions for each iteration
# Param: features - a numerical matrix of features
# Param: outcomes - a numerical vector of outcomes
# Param: categorical - a vector representing which features are categorical
# Param: iterations to call appropriate functions
# Return: classifier consisting of a linear combination of decision stumps
make_classifier <- function(features, outcomes, categorical, iterations, verbose) {

  # PREPARE INPUT
  # --------------------------------------------------------------------------------
  ordered_index <- matrix(NA, nrow = nrow(features), ncol = ncol(features))
  for (i in 1:ncol(features)) {
    ordered_index[, i] <- order(features[, i]) - 1
  }

  # CALL C++ CODE
  # --------------------------------------------------------------------------------
  classifier <- adaboost(features, ordered_index, outcomes, categorical, iterations, verbose)

  return(classifier)
}


