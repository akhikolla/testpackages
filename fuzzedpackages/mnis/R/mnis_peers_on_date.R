
#' A tibble with information on all peers who were members of the House of Commons on the date specificed (if only date1 is included as a parameter), or on or between the two dates if both date1 and date2 are specified.
#' @param date1 The date to return the list of peers from. Defaults to current system date. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}.
#' @param date2 An optional query parameter. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. If not NULL, the function returns a list of all peers in the House of Lords between date2 and date1. Defaults to NULL.
#' @param tidy Fix the variable names in the tibble to remove extra characters, superfluous text and convert variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#'
#' @return A tibble with information on all Peers who were members of the House of Lords on the date specificed (if only date1 is included as a parameter), or on or between the two dates if both date1 and date2 are specified.
#' @export
#' @seealso \code{\link{mnis_party_state}} \code{\link{mnis_peers_on_date}}
#' @examples \dontrun{
#'
#' x <- mnis_peers_on_date("2017-01-01")
#'
#' }

mnis_peers_on_date <- function(date1 = Sys.Date(), date2=NULL, tidy = TRUE, tidy_style="snake_case"){

  baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/House=Lords|Membership=all|lordsmemberbetween="

  date1 <- as.Date(date1)

  if(is.null(date2)==FALSE) {
    date2 <- as.Date(date2)
  }

  if(is.null(date2)==TRUE) {
    date2 <- date1
  } else if (date1 > date2) {
    date3 <- date1
    date1 <- date2
    date2 <- date3
    rm(date3)
  }

  query <- paste0(baseurl,date1,"and",date2,"/")

  got <- httr::GET(query, httr::accept_json())

  if (httr::http_type(got) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }

  got <- mnis::tidy_bom(got)

  got <- jsonlite::fromJSON(got, flatten = TRUE)

  lords <- got$Members$Member

  lords <- tibble::as_tibble(lords)

  if (tidy == TRUE) {

    lords <- mnis::mnis_tidy(lords, tidy_style)

    lords

  } else {

    lords

  }

}
