#' Temperature datasets from Canada.
#'
#' Temperature time series from 35 weather stations from Canada.  This dataset
#' is a classic one and was used in famous package \code{fda}. We have made a few
#' changes in this dataset.
#'
#' @usage
#' data("datasetCanada")
#'
#' @format A list with two entries: \code{m_cood} and \code{m_data}.
#' \describe{
#'   \item{\code{m_coord}}{a \code{tibble} with latitude, logitude and the name
#'   of stations.}
#'   \item{m_data}{a \code{tibble} where each column is the time series from a
#'   weather station.}
#' }
#' @source the \code{CanadianWeather} dataset from the \code{R} package
#' \code{fda}.
"datasetCanada"
