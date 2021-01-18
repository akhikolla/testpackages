
#' Calls the API to return a tibble with details on the number of Lords and their affiliations.
#' @param date Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. The API will return data on the state of the House of Lords on that date. Defaults to the current system date.
#' @param tidy If TRUE, fixes the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @return A tibble with information on the numbers of different types of Lords on a given date.
#' @keywords mnis
#' @export
#' @seealso \code{\link{mnis_reference}}
#' @examples \dontrun{
#'
#' x <- mnis_lords_type()
#'
#' }

mnis_lords_type <- function(date = Sys.Date(), tidy = TRUE, tidy_style="snake_case") {

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/LordsByType/"

    date <- as.Date(date)

    query <- paste0(baseurl, date, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- tibble::as_tibble(as.data.frame(got$LordsByType))

    if (tidy == TRUE) {

        x <- mnis::mnis_tidy(x, tidy_style)

        x

    } else {

        x

    }

}

