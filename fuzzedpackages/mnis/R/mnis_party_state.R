

#' A tibble with information on the numbers and gender of MPs, by party, for the given date.
#' @param house The house of parliament. Accepts either 'Lords' or 'Commons'. Defaults to 'Commons'.
#' @param date Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to the current system date.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @return A tibble with information on the numbers and gender of MPs, by party, by party, for the given date.
#' @keywords mnis
#' @seealso \code{\link{mnis_mps_on_date}} \code{\link{mnis_peers_on_date}}
#' @export
#' @examples \dontrun{
#'
#' x <- mnis_party_state('2012-01-12')
#'
#' }

mnis_party_state <- function(house = "Commons", date = Sys.Date(), tidy = TRUE, tidy_style="snake_case") {

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/houseOverview/"

    query <- paste0(baseurl, house, "/", date, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- tibble::as_tibble(got$HouseOverview$Party)

    if (tidy == TRUE) {

        x <- mnis_tidy(x, tidy_style)

    } else {

        x

    }

}
