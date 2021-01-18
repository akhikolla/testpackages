#' pandemicPredicted objects: Predictions made from a fitted PandemicLP model
#'
#' The \pkg{PandemicLP} prediction function returns an object of S3 class
#' \code{pandemicPredicted}, which is a list containing the components described below. \cr
#' \cr
#'
#' @name pandemicPredicted-objects
#'
#' @section Elements for \code{pandemicPredicted} objects:
#' \describe{
#'   \item{\code{predictive_Long}}{
#'   The full sample of the predictive distribution for the long-term prediction.
#'   The prediction is for daily new cases.
#'   }
#'   \item{\code{predictive_Short}}{
#'   The full sample of the predictive distribution for the short-term prediction.
#'   The prediction is for daily cumulative cases.
#'   }
#'   \item{\code{data}}{
#'   The data passed on from the \code{\link{pandemicEstimated-objects}} under the element \code{Y$data}.
#'   }
#'   \item{\code{location}}{
#'   A string with the name of the location.
#'   }
#'   \item{\code{cases_type}}{
#'   A string with either "confirmed" or "deaths" to represent the type of data that has been fitted and predicted.
#'   }
#'   \item{\code{pastMu}}{
#'   The fitted means of the data for the observed data points.
#'   }
#'   \item{\code{futMu}}{
#'   The predicted means of the data for the predicted data points.
#'   }
#' }
#'
NULL
