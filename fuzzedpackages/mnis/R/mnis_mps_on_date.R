
#' A tibble with information on all MPs who were members of the House of Commons on the date specificed (if only date1 is included as a parameter), or on or between the two dates if both date1 and date2 are specified.
#' @param date1 The date to return the list of mps from. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to current system date.
#' @param date2 An optional query parameter. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. If not NULL, the function returns a list of all MPs who were members between date2 and date1. Defaults to NULL.
#' @param tidy Fix the variable names in the tibble to remove extra characters, superfluous text and convert variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#'
#' @return A tibble with information on all MPs who were members of the House of Commons on the date specificed (if only date1 is included as a parameter), or on or between the two dates if both date1 and date2 are specified.
#' @export
#'
#' @examples \dontrun{
#'
#' x <- mnis_mps_on_date()
#'
#' }

mnis_mps_on_date <- function(date1 = Sys.Date(), date2=NULL, tidy = TRUE, tidy_style="snake_case"){

  baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/House=Commons|Membership=all|commonsmemberbetween="

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

  mps <- got$Members$Member

  mps <- tibble::as_tibble(mps)

  if(.Platform$OS.type=="windows"){

    mps$MemberFrom <- stringi::stri_trans_general(mps$MemberFrom, "latin-ascii")

    mps$MemberFrom <- gsub("Ynys MA\U00B4n", "Ynys M\U00F4n", mps$MemberFrom)

  }

  if (tidy == TRUE) {

    mps <- mnis::mnis_tidy(mps, tidy_style)

    mps

  } else {

    mps

  }

}
