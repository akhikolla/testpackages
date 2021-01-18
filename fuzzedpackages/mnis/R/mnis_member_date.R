
#' Returns a tibble with a member's status on a given date.
#' @param ID The ID of the member. Currently this only accepts IDs from the default membership ID scheme. If empty, the function stops and no data is returned.
#' @param date Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Return details on the requested member's status on that date. Defaults to the current system date.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case"
#' @return Returns a tibble with the given member's status on the given date.
#' @keywords mnis
#' @export
#' @seealso \code{\link{mnis_mps_on_date}}
#' @examples \dontrun{
#'
#' x <- mnis_member_date(172)
#'
#' }

mnis_member_date <- function(ID = NULL, date = Sys.Date(), tidy = TRUE, tidy_style="snake_case") {

    if (missing(ID)) {
        stop("The ID parameter cannot be empty, please specify a Member of Parliament or a Peer.")
    }

    ID <- as.character(ID)

    date <- as.Date(date)

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/member/historical/"

    query <- paste0(baseurl, ID, "/", date, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- as.list(got$Member)

    x <- unlist(x)

    x <- t(x)

    x <- tibble::as_tibble(x)

    if (tidy == TRUE) {

        x <- mnis::mnis_tidy(x, tidy_style)

        x

    } else {

        x

    }

}
