#' pandemicEstimated objects: Fitted PandemicLP model
#'
#' The \pkg{PandemicLP} model-fitting functions return an object of S3 class
#' \code{pandemicEstimated}, which is a list containing the components described below.
#'
#'
#' @name pandemicEstimated-objects
#'
#' @section Elements for \code{pandemicEstimated} objects:
#' \describe{
#'   \item{\code{model_name}}{
#'   The model name used.
#'   }
#'   \item{\code{multiwaves}}{
#'   number of waves.
#'   }
#'   \item{\code{seasonal_effect}}{
#'   String vector of the  weekdays' name with sazonal effect. Only for one wave model.
#'   }
#'   \item{\code{cases.type}}{
#'   The type of cases of interest used in modeling the epidemic.
#'   }
#'   \item{\code{config.inputs}}{
#'   A list with the main input arguments used in \code{\link{pandemic_model}}.
#'   }
#'   \item{\code{priors}}{
#'   A list with information on the prior distributions used and model restrictions (if there are any).
#'   }
#'    \item{\code{fit}}{
#'   An object of S4 Class \code{\link[rstan]{stanfit}} representing the fitted results via
#'   \code{\link[rstan]{sampling}}. For additional information about this element \code{fit},
#'   see \code{\link[rstan]{stanfit}}.
#'   }
#'   \item{\code{Y}}{
#'   A list with the \code{data}.
#'   }
#' }
#'
NULL
