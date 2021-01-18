#' mnis: Easy Downloading Capabilities for the Member Name Information Service
#'
#' An API package for the Members' Name Information Service operated by the UK parliament. The package is intended to simplify pulling data from an API for users unfamiliar with APIs. Documentation for the API itself can be found here: http://data.parliament.uk/membersdataplatform/default.aspx.
#'
#' The package includes a built in function to remove a byte-order mark from the API data, and a parameter \code{tidy} with each function that converts variable names into an R friendly style, removing non-alphanumeric characters and converting to snake_case when equal to TRUE, its default value.
#'
#' All functions requests data in JSON format and parse it to a tibble. The exception is \code{\link{mnis_constituency_results}} which returns a single object containing a list (with constituency details) and a tibble (with election results).
#'
#' None of the functions included in this package appear to run into the API rate limit, although there may be restrictions to custom requests, which allow a maximum of three parameters.
#'
#'
#' @section mnis functions:
#'
#' \code{\link{mnis_additional}}
#'
#' \code{\link{mnis_all_members}}
#'
#' \code{\link{mnis_base}}
#'
#' \code{\link{mnis_constituency_results}}
#'
#' \code{\link{mnis_department}}
#'
#' \code{\link{mnis_eligible}}
#'
#' \code{\link{mnis_extra}}
#'
#' \code{\link{mnis_full_biog}}
#'
#' \code{\link{mnis_general_election_results}}
#'
#' \code{\link{mnis_joined_between}}
#'
#' \code{\link{mnis_lords_type}}
#'
#' \code{\link{mnis_member_date}}
#'
#' \code{\link{mnis_party_state}}
#'
#' \code{\link{mnis_reference}}
#'
#' @docType package
#' @name mnis
#' @import utils
#' @import httr
#' @import jsonlite
#' @import dplyr
#' @import tibble
#' @import stringi
#' @useDynLib mnis, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
